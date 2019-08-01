// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
#include "Doussan.h";

#include "Organism.h"
#include "soil.h"



namespace CRootBox {



void Doussan::init(int ot)
{
	double simtime = plant.getSimTime();
	nodes = plant.getNodes();
	segs =plant.getSegments(ot);
	auto segsO =plant.getSegmentOrigins(ot);
	auto segsCT =plant.getSegmentCTs(ot);
	kx_.resize(segs.size());
	kr_.resize(segs.size());
	for (int i=0; i<segs.size(); i++) {
		int subType = segsO[i]->getParam()->supType;
		double age = simtime - segsCT[i];
		kx_[i] = axialConductivity(ot, subType, age);
		kr_[i] = radialConductivity(ot, subType, age);
		auto n1 = nodes[segs[i].x];
		auto n2 = nodes[segs[i].y];
// 		auto v = n2.minus(n1); v2, length
	}
}

std::vector<double> Doussan::getI()
{
	std::vector<double> I = std::vector<double>(4*segs.size());
	for (int i=0; i<segs.size(); i+=4) {
		int i = segs[i].x;
		int j = segs[i].y;
		I[i] = i;
		I[i+1] = i;
		I[i+2] = j;
		I[i+3] = j;
	}
	return I;
}

std::vector<double> Doussan::getJ()
{
	std::vector<double> J = std::vector<double>(4*segs.size());
	for (int i=0; i<segs.size(); i+=4) {
		int i = segs[i].x;
		int j = segs[i].y;
		J[i] = i;
		J[i+1] = j;
		J[i+2] = j;
		J[i+3] = i;
	}
	return J;
}

std::vector<double> Doussan::getV()
{
	std::vector<double> V = std::vector<double>(4*segs.size());
	for (int i=0; i<segs.size(); i+=4) {

	}
	return V;
}

std::vector<double> Doussan::getB()
{
	std::vector<double> b = std::vector<double>(4*segs.size());
	for (int i=0; i<segs.size(); i+=4) {
		int i = segs[i].x;
		int j = segs[i].y;
		auto n1 = nodes[i];
		auto n2 = nodes[j];
		auto v = n2.minus(n1).times(1./length[i]);


	}
	return b;
}


} // namespace
