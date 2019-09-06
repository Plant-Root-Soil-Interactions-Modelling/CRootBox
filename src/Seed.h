#ifndef SEED_H_
#define SEED_H_

#include "Organ.h"
#include "Organism.h"
#include "seedparameter.h"

namespace CRootBox {

/**
 * Seed
 *
 * simulate calls the simulate method of the stem, and base roots
 */
class Seed : public Organ
{
public:

	Seed(Organism* plant) :Organ(plant, nullptr, Organism::ot_seed, 0, 0.) { };
	virtual ~Seed() { }; // TODO delete all

	virtual int organType() const override { return Organism::ot_seed; }

	void initialize();

	void simulate(double dt, bool verbose = false) override; // TODO

    SeedSpecificParameter* param() const { return (SeedSpecificParameter*)param_; }

	virtual std::string toString() const override;

	int getNumberOfRootCrowns() const { return numberOfRootCrowns; }

	// default positions
	int basalType = 4;
	int shootborneType = 5;
	int tillerType = 4;

protected:

    int numberOfRootCrowns = 0;
	int getParamSubType(int organtype, std::string str);


};

}

#endif /* Seed_H_ */
