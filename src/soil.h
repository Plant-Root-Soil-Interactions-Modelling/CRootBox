// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
#ifndef SOIL_H
#define SOIL_H

#include "mymath.h"
#include "sdf.h"
#include <cmath>
#include <limits>


namespace CRootBox  {

    class Root;

    /**
     * Base class to look up for a scalar soil property
     */
    class SoilLookUp
    {
    public:

        const double inf = std::numeric_limits<double>::infinity();

        SoilLookUp() { };
        virtual ~SoilLookUp() { };

        /**
         * Returns a scalar property of the soil scaled from 0..1.
         * if you want include periodicity, use periodic(pos) instead of pos
         *
         * @param pos       position [cm], (normally, root->getNode(root->getNumberOfNodes()-1))
         * @param root      the root that wants to know the scalar property
         *                  in some situation this might be usefull (e.g. could increase look up speed from a unstructured mesh)
         * \return          scalar soil property
         */
        virtual double getValue(const Vector3d& pos, const Root* root = nullptr) const { return 1.; } ///< Returns a scalar property of the soil, 1. per default

        virtual std::string toString() const { return "SoilLookUp base class"; } ///< Quick info about the object for debugging

        /**
         * sets the periodic boundaries, periodicity is used if bounds are set, and if it is supportet by
         * the specialised class
         */
        void setPeriodicDomain(double minx_, double maxx_, double miny_ , double maxy_, double minz_, double maxz_) {
            periodic_ = true;
            minx = minx_;
            xx = maxx_-minx;
            miny = miny_;
            yy = maxy_-miny;
            minz = minz_;
            zz = maxz_-minz;
        }

        void setPeriodicDomain(double minx, double maxx, double miny , double maxy) {
            this->setPeriodicDomain(minx, maxx, miny, maxy, 0, inf );
        }

        void setPeriodicDomain(double minx, double maxx) {
            this->setPeriodicDomain(minx, maxx, -inf, inf, 0, inf );
        }

        Vector3d periodic(const Vector3d& pos) {  //< maps point into periodic domain
            if (periodic_) {
                    Vector3d p = pos;
                    if (!std::isinf(xx)) { // periodic in x
                            p.x -= minx;
                            p.x = (p.x/xx - (int)(p.x/xx))*xx;
                            p.x += minx;
                    }
                    if (!std::isinf(yy)) { // periodic in y
                            p.y -= miny;
                            p.y = (p.y/yy - (int)(p.y/yy))*yy;
                            p.y += miny;
                    }
                    if (!std::isinf(zz)) { // periodic in z
                            p.z -= minz;
                            p.z = (p.z/zz - (int)(p.z/zz))*zz;
                            p.z += minz;
                    }
                    return p;
            } else {
                    return pos;
            }
        }

    private:
        bool periodic_ = false;
        double minx=0., xx=0., miny=0., yy=0, minz=0., zz=0.;

    };



    /**
     * Looks up a value based on a signed distance function
     */
    class SoilLookUpSDF : public SoilLookUp
    {
    public:
        SoilLookUpSDF(): SoilLookUpSDF(nullptr) { } ///< Default constructor

        /**
         * Creaets the soil property from a signed distance function,
         * inside the geometry the value is largest
         *
         * @param sdf_      the signed distance function representing the geometry
         * @param max_      the maximal value of the soil property
         * @param min_      the minimal value of the soil property
         * @param slope_    scales the linear gradient of the sdf (note that |grad(sdf)|= 1)
         */
        SoilLookUpSDF(SignedDistanceFunction* sdf_, double max_=1, double min_=0, double slope_=1) {
            this->sdf=sdf_;
            fmax = max_;
            fmin = min_;
            slope = slope_;
        } ///< Creates the soil property from a signed distance function

        virtual ~SoilLookUpSDF() { };

        virtual double getValue(const Vector3d& pos, const Root* root = nullptr) const override {
            double c = -sdf->getDist(pos)/slope*2.; ///< *(-1), because inside the geometry the value is largest
            c += (fmax-fmin)/2.; // thats the value at the boundary
            return std::max(std::min(c,fmax),fmin);
        } ///< returns fmin outside of the domain and fmax inside, and a linear ascend according slope

        virtual std::string toString() const override { return "SoilLookUpSDF"; } ///< Quick info about the object for debugging

        SignedDistanceFunction* sdf; ///< signed distance function representing the geometry
        double fmax; ///< maximum is reached within the geometry at the distance slope
        double fmin; ///< minimum is reached outside of the geometry at the distance slope
        double slope; ///< half length of linear interpolation between fmax and fmin
    };



    /**
     * Scales the root elongation with fixed value same for each root
     */
    class ProportionalElongation : public SoilLookUp
    {
    public:

        void setScale(double s) { scale = s; }

        void setBaseLookUp(SoilLookUp* baseLookUp) { this->baseLookUp=baseLookUp; } ///< proportionally scales a base soil look up

        virtual double getValue(const Vector3d& pos, const Root* root = nullptr) const override {
            if (baseLookUp==nullptr) {
                    return scale;
            } else {
                    return baseLookUp->getValue(pos,root)*scale;  //super impose scaling on a base soil look up function
            }
        }

        virtual std::string toString() const override { return "ProportionalElongation"; } ///< Quick info about the object for debugging

    protected:
        double scale = 1.;
        SoilLookUp* baseLookUp = nullptr;

    };



    /**
     * 1D look up table
     */
    class Grid1D  : public SoilLookUp
    {
    public:

        Grid1D() {
            n=0;
            data = std::vector<double>(0);
            grid = std::vector<double>(1);
        }

        Grid1D(size_t n, std::vector<double>grid, std::vector<double> data): n(n), grid(grid), data(data) {
            assert(grid.size()==n);
            assert(data.size()==n-1);
        };

        virtual size_t map(double x) const {
            unsigned int jr,jm,jl;
            jl = 0;
            jr = n - 1;
            while (jr-jl > 1) {
                    jm=(jr+jl) >> 1; // thats a divided by two
                    if (x >= grid[jm])
                        jl=jm;
                    else
                        jr=jm;
            }
            return jl;
        } ///< Generic way to perform look up in an ordered table, overwrite by faster method if appropriate

        virtual double getValue(const Vector3d& pos, const Root* root = nullptr) const override {
            size_t i = map(pos.z);
            i = std::max(int(i),0);
            i = std::min(int(i),int(n-1));
            return data[i];
        } ///< Returns the data of the 1d table, repeats first or last entry if out of bound

        virtual std::string toString() const override { return "RectilinearGrid1D"; } ///< Quick info about the object for debugging

        size_t n;
        std::vector<double> grid;
        std::vector<double> data;
    };



    /**
     *  1D look up table with equidistant spacing
     */
    class EquidistantGrid1D : public Grid1D
    {
    public:

        EquidistantGrid1D(double a, double b, size_t n): a(a), b(b) {
            this->n =n;
            makeGrid(a,b,n);
            this->data = std::vector<double>(n-1);
        }

        EquidistantGrid1D(double a, double b, const std::vector<double>& data): a(a), b(b) {
            this->n = data.size()+1;
            makeGrid(a,b,n);
            this->data = data;
        }

        void makeGrid(double a, double b, size_t n) {
            this->grid = std::vector<double>(n);
            for (size_t i=0; i<n; i++) {
                    grid[i] = a + (b-a)/double(n-1)*i;
            }
        }

        virtual size_t map(double x) const override {
            return std::floor((x-a)/(b-a)*(n-1));
        }

        virtual std::string toString() const  override{ return "LinearGrid1D"; } ///< Quick info about the object for debugging

        double a;
        double b;

    };


    /**
     * Fantastic
     */
    class RectilinearGrid3D  : public SoilLookUp
    {
    public:

        RectilinearGrid3D(Grid1D* xgrid, Grid1D* ygrid, Grid1D* zgrid) :xgrid(xgrid), ygrid(ygrid), zgrid(zgrid) {
            nx = xgrid->n;
            ny = ygrid->n;
            nz = zgrid->n;
            data = std::vector<double>(nx*ny*nz);
        };

        virtual ~RectilinearGrid3D() { };

        virtual size_t map(double x, double y, double z) const {
            size_t i = xgrid->map(x);
            size_t j = ygrid->map(y);
            size_t k = zgrid->map(z);
            return i*(nx*ny)+j*ny+k; // or whatever
        }

        double getData(size_t i, size_t j, size_t k) {
            return data.at(map(i,j,k));
        }

        void setData(size_t i, size_t j, size_t k, double d) {
            data.at(map(i,j,k)) = d;
        }

        Grid1D* xgrid;
        Grid1D* ygrid;
        Grid1D* zgrid;

        size_t nx,ny,nz;
        std::vector<double> data;

    };

    /**
     *  Even more fantastic
     */
    class EquidistantGrid3D : public RectilinearGrid3D
    {
        EquidistantGrid3D(double x0, double xe, int nx, double y0, double ye, int ny, double z0, double ze, int nz) : RectilinearGrid3D(new EquidistantGrid1D(x0,xe,nx),new EquidistantGrid1D(y0,ye,ny),new EquidistantGrid1D(z0,ze,nz)) {
        }

        virtual ~EquidistantGrid3D() {
            delete xgrid;
            delete ygrid;
            delete zgrid;
        }

    };

} // end namespace CRootBox

#endif
