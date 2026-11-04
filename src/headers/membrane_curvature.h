#ifndef MOSAICS_MEMBRANE_CURVATURE_H
#define MOSAICS_MEMBRANE_CURVATURE_H
// Numerical kernel only: no trajectory storage, MPI, or MOSAICS global state.
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>
namespace membrane_curvature
{
    struct Point
    {
        double x, y, z;
    };
    struct Curvature
    {
        double mean, gaussian;
        bool valid;

        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //                                                                                                           //
        // This constructor initializes an unsupported curvature sample.                                             //
        //                                                                                                           //
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
        Curvature() : mean(0), gaussian(0), valid(false)
        {
        }
    };

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // This function places a coordinate inside its periodic box interval.                                       //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    inline double wrap(double x, double length)
    {
        return x - length * std::floor(x / length);
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // This function returns the minimum-image separation along one box dimension.                               //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    inline double minimum_image(double x, double length)
    {
        return x - length * std::floor(x / length + 0.5);
    }

    // n = (-h_x,-h_y,1)/sqrt(1+h_x^2+h_y^2). Use the convention in which
    // a dome viewed from the positive normal has positive mean curvature.
    // normal_sign = +1: common +Z normal; -1: reverse the normal (H only).


    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // This function computes mean and Gaussian curvature from the surface derivatives.                          //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    inline Curvature from_derivatives(double hx, double hy, double hxx,
        double hxy, double hyy, double normal_sign)
    {
        Curvature c;
        const double q = 1.0 + hx*hx + hy*hy;
        c.mean = -normal_sign * ((1.0+hy*hy)*hxx - 2.0*hx*hy*hxy
            + (1.0+hx*hx)*hyy) / (2.0*q*std::sqrt(q));
        c.gaussian = (hxx*hyy - hxy*hxy)/(q*q);
        c.valid = std::isfinite(c.mean) && std::isfinite(c.gaussian);
        return c;
    }
    class Surface
    {
        public:
            int nx, ny;
            double lx, ly, dx, dy;
            bool periodic;
            std::vector<double> height;
            std::vector<unsigned char> supported;
            std::vector<Curvature> curvature;

            ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
            //                                                                                                           //
            // This constructor allocates a surface grid with the requested geometry and spacing.                        //
            //                                                                                                           //
            ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
            Surface(double length_x, double length_y, double spacing, bool pbc)
            : lx(length_x), ly(length_y), periodic(pbc)
            {
                if (!(lx > 0 && ly > 0 && spacing > 0) ||
                    !std::isfinite(lx+ly+spacing))
                {
                    throw std::runtime_error("Invalid surface dimensions or spacing.");
                }
                const double count_x = std::ceil(lx/spacing);
                const double count_y = std::ceil(ly/spacing);
                // Explicit bound prevents integer overflow and accidental huge allocations.
                if (count_x > 100000 || count_y > 100000 || count_x*count_y > 20000000)
                {
                    throw std::runtime_error("Surface exceeds 20 million nodes; increase -surf_step.");
                }
                nx = static_cast<int>(count_x) + (periodic ? 0 : 1);
                ny = static_cast<int>(count_y) + (periodic ? 0 : 1);
                if (nx < 5 || ny < 5)
                {
                    throw std::runtime_error("Surface needs at least five nodes in each direction.");
                }
                dx = lx/(periodic ? nx : nx-1);
                dy = ly/(periodic ? ny : ny-1);
                const std::size_t size = static_cast<std::size_t>(nx)*ny;
                height.resize(size);
                supported.resize(size);
                curvature.resize(size);
            }

            ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
            //                                                                                                           //
            // This function returns the surface array index, applying periodic wrapping when requested.                 //
            //                                                                                                           //
            ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
            std::size_t at(int x, int y) const
            {
                if (periodic)
                {
                    x = (x%nx+nx)%nx;
                    y = (y%ny+ny)%ny;
                }
                return static_cast<std::size_t>(y)*nx+x;
            }

            // Gaussian kernel regression of Z on XY, using one surface atom per lipid.
            // The finite Gaussian support is a square of half-width radius (a separable
            // kernel), not a circular stamp. Use >=4 sigma (the default support is 4 sigma).
            // The normalization prevents an empty grid cell from being treated as z=0.
            // Neff=(sum w)^2/sum(w^2) and sum w mask weakly supported regions.


            ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
            //                                                                                                           //
            // This function computes a Gaussian-smoothed height surface and its support mask.                           //
            //                                                                                                           //
            ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
            void reconstruct(const std::vector<Point>& points, double sigma, double radius,
                double min_weight, double min_neff, double lz)
            {
                if (!(sigma > 0 && radius > 0 && min_weight > 0 && min_neff >= 1 && lz > 0))
                {
                    throw std::runtime_error("Invalid smoothing/support parameters or Z box.");
                }
                if (periodic && radius >= 0.5*std::min(lx, ly))
                {
                    throw std::runtime_error("-smooth_r must be less than half the shortest XY box length.");
                }
                std::fill(height.begin(), height.end(), 0.0);
                std::fill(supported.begin(), supported.end(), 0);
                std::fill(curvature.begin(), curvature.end(), Curvature());
                if (points.empty())
                {
                    return;
                }
                std::vector<double> sum_w(height.size(), 0), sum_w2(height.size(), 0);
                const double pi = 3.14159265358979323846;
                double sum_sin = 0, sum_cos = 0;
                for (const Point& p : points)
                {
                    if (!std::isfinite(p.x+p.y+p.z))
                    {
                        throw std::runtime_error("Nonfinite surface coordinate.");
                    }
                    sum_sin += std::sin(2*pi*wrap(p.z, lz)/lz);
                    sum_cos += std::cos(2*pi*wrap(p.z, lz)/lz);
                }

                // A common image per leaflet avoids artificial steps when Z wraps.
                if (std::hypot(sum_sin, sum_cos) < 0.1*points.size())
                {
                    throw std::runtime_error("Leaflet Z positions do not define a single thin height surface.");
                }
                const double anchor = wrap(std::atan2(sum_sin, sum_cos)*lz/(2*pi), lz);
                double z_min = std::numeric_limits<double>::infinity(), z_max = -z_min;
                for (const Point& p : points)
                {
                    const double z = anchor+minimum_image(p.z-anchor, lz);
                    z_min = std::min(z_min, z);
                    z_max = std::max(z_max, z);
                }
                if (z_max-z_min >= 0.5*lz)
                {
                    throw std::runtime_error("Leaflet Z span "+std::to_string(z_max-z_min)+
                        " nm exceeds half its Z box. Check reference leaflet assignments, "
                        "flip-flop and unwrapping; a single height surface is required.");
                }
                for (const Point& p : points)
                {
                    const double x = periodic ? wrap(p.x, lx) : p.x;
                    const double y = periodic ? wrap(p.y, ly) : p.y;
                    const double z = anchor+minimum_image(p.z-anchor, lz);
                    if(!periodic && (x < -radius || x > lx+radius || y < -radius || y > ly+radius))
                    {
                        continue;
                    }
                    int x0 = static_cast<int>(std::ceil((x-radius)/dx));
                    int x1 = static_cast<int>(std::floor((x+radius)/dx));
                    int y0 = static_cast<int>(std::ceil((y-radius)/dy));
                    int y1 = static_cast<int>(std::floor((y+radius)/dy));
                    if (!periodic)
                    {
                        x0 = std::max(0, x0);
                        x1 = std::min(nx-1, x1);
                        y0 = std::max(0, y0);
                        y1 = std::min(ny-1, y1);
                    }
                    for (int j = y0;j <= y1;++j)
                    {
                        for(int i = x0;i <= x1;++i)
                        {
                            const double rx = i*dx-x, ry = j*dy-y, r2 = rx*rx+ry*ry;
                            // Separable Gaussian support: |rx|, |ry| <= radius.
                            const double w = std::exp(-0.5*r2/(sigma*sigma));
                            const std::size_t k = at(i, j);
                            height[k] += w*(z-anchor);
                            sum_w[k] += w;
                            sum_w2[k] += w*w;
                        }
                    }
                }
                for (std::size_t k = 0;k<height.size();++k)
                {
                    if (sum_w[k] >= min_weight && sum_w[k]*sum_w[k] >= min_neff*sum_w2[k])
                    {
                        height[k] = height[k]/sum_w[k]+anchor;
                        supported[k] = 1;
                    }
                }
            }

            // Second-order central differences. Require the whole 3x3 stencil to
            // be supported, including diagonals for h_xy; never bridge invalid holes.


            ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
            //                                                                                                           //
            // This function computes curvature from supported first and second surface derivatives.                     //
            //                                                                                                           //
            ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
            void differentiate(double normal_sign)
            {
                std::fill(curvature.begin(), curvature.end(), Curvature());
                for (int j = 0;j<ny;++j)
                {
                    for(int i = 0;i<nx;++i)
                    {
                        if(!periodic && (i == 0 || j == 0 || i == nx-1 || j == ny-1))
                        {
                            continue;
                        }
                        bool valid = true;
                        for(int b = -1;b <= 1;++b)
                        {
                            for(int a = -1;a <= 1;++a)
                            {
                                if (!supported[at(i+a, j+b)])
                                {
                                    valid = false;
                                }
                            }
                        }
                        if (!valid)
                        {
                            continue;
                        }
                        const double h = height[at(i, j)];
                        const double hx = (height[at(i+1, j)]-height[at(i-1, j)])/(2*dx);
                        const double hy = (height[at(i, j+1)]-height[at(i, j-1)])/(2*dy);
                        const double hxx = (height[at(i+1, j)]-2*h+height[at(i-1, j)])/(dx*dx);
                        const double hyy = (height[at(i, j+1)]-2*h+height[at(i, j-1)])/(dy*dy);
                        const double hxy = (height[at(i+1, j+1)]-height[at(i+1, j-1)]
                            -height[at(i-1, j+1)]+height[at(i-1, j-1)])/(4*dx*dy);
                        curvature[at(i, j)] = from_derivatives(hx, hy, hxx, hxy, hyy, normal_sign);
                    }
                }
            }

            ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
            //                                                                                                           //
            // This function interpolates both curvatures at the requested mapping-atom XY position.                     //
            //                                                                                                           //
            ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
            Curvature sample(double x, double y) const
            {
                Curvature result;
                if (!std::isfinite(x+y))
                {
                    return result;
                }
                if(periodic)
                {
                    x = wrap(x, lx);
                    y = wrap(y, ly);
                }
                else if(x<0 || x >= lx || y<0 || y >= ly)
                {
                    return result;
                }
                const double gx = x/dx, gy = y/dy;
                const int i = static_cast<int>(std::floor(gx)), j = static_cast<int>(std::floor(gy));
                const double a = gx-i, b = gy-j;
                result.valid = true;
                for(int q = 0;q<2;++q)
                {
                    for(int p = 0;p<2;++p)
                    {
                        const double w = (p?a:1-a)*(q?b:1-b);
                        if(w <= 1e-14)
                        {
                            continue;
                        }
                        const Curvature& c = curvature[at(i+p, j+q)];
                        if(!c.valid)
                        {
                            return Curvature();
                        }
                        result.mean += w*c.mean;
                        result.gaussian += w*c.gaussian;
                    }
                }
                return result;
            }
    };
}

// namespace membrane_curvature
#endif
