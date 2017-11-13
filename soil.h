#ifndef SOIL_H
#define SOIL_H

#include "mymath.h"
#include <cmath>
#include "sdf.h"

class Root;



/**
 * Look up for a scalar soil property
 */
class SoilLookUp
{
public:
	SoilLookUp() { };
	virtual ~SoilLookUp() { };

	/**
	 * Returns a scalar property of the soil scaled from 0..1
	 *
	 * @param pos       position [cm], (normally, root->getNode(root->getNumberOfNodes()-1))
	 * @param root      the root that wants to know the scalar property
	 *                  in some situation this might be usefull (e.g. could increase look up speed from a unstructured mesh)
	 * \return          scalar soil property
	 */
	virtual double getValue(const Vector3d& pos, const Root* root = nullptr) const { return 1.; } ///< Returns a scalar property of the soil, 1. per default

	virtual std::string toString() const { return "SoilLookUp base class"; } ///< Quick info about the object for debugging

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

	virtual double getValue(const Vector3d& pos, const Root* root = nullptr) const override {
		return scale;
	}

	virtual std::string toString() const override { return "ProportionalElongation"; } ///< Quick info about the object for debugging

protected:
	double scale = 1.;

};



/**
 * 1D look up table
 *
 *
 */
class RectilinearGrid1D  : public SoilLookUp
{
public:

	RectilinearGrid1D() {
		n=0;
		data = std::vector<double>(0);
		grid = std::vector<double>(1);
	}

	RectilinearGrid1D(size_t n, std::vector<double>grid, std::vector<double> data): n(n), grid(grid), data(data) {
		assert(grid.size()==n);
		assert(data.size()==n-1);
	};

	virtual size_t map(double x) const {
		unsigned int jr,jm,jl;
		jl=0;
		jr=n;
		while (jr-jl > 1) {
			jm=(jr+jl) >> 1; // thats divided by two
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
class EquidistantGrid1D : public RectilinearGrid1D
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

	RectilinearGrid3D(RectilinearGrid1D* xgrid, RectilinearGrid1D* ygrid, RectilinearGrid1D* zgrid) :xgrid(xgrid), ygrid(ygrid), zgrid(zgrid) {
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

	RectilinearGrid1D* xgrid;
	RectilinearGrid1D* ygrid;
	RectilinearGrid1D* zgrid;

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
	};

};



#endif
