#include "ClusterLabelling.h"


ClusterSystem::Particle::Particle(const int id, const double theta) :
	_id(id), _theta(theta)
{}

ClusterSystem::Contact::Contact(Particle* p1, Particle* p2, 
	double fn, double fs, double f) :
		_particle1(p1), _particle2(p2),
	_normalForce(fn), _shearForce(fs), _totalForce(f)
{}

void ClusterSystem::Contact::calculateRelativeAngle()
{
	_dTheta = abs(_particle1->_theta - _particle2->_theta);
	_dTheta = fmod(_dTheta, M_PI);
	if (_dTheta > M_PI_2) _dTheta = M_PI - _dTheta;
}


ClusterSystem::ClusterAnalysis::ClusterAnalysis(const vector<BaseParticle*>& P, const double accuracy)
{
	// Generate particles from DEM particles...
	_particles.reserve(P.size());
	for (const auto& p : P) _particles.push_back(Particle(p->getID(), p->getRot()));

	// Generate contacts from DEM contacts...
	_contacts.reserve(P.size() * 10);
	for (const auto& p : P) {
		int p1 = p->getID();
		for (const auto& c0 : p->getContactInfo()) {
			int p2 = c0.first;
			const auto& c = c0.second;
			_contacts.push_back(Contact(
				&_particles[p1], 
				&_particles[p2],
				c._normalForce.getLength(),
				c._shearForce.getLength(),
				c._force.getLength()));
		}
	}

	// calculate threshold vector
	double dtheta = 0.0;
	while (dtheta < M_PI_2) {
		_dTheta0.push_back(dtheta);
		dtheta += accuracy;
	}
	_mclus.resize(_dTheta0.size(), vector<int>(_particles.size(), -1));
	_nclus.resize(_dTheta0.size(), vector<int>(_particles.size(), -1));
	_correlation.resize(_dTheta0.size(), 0.0);
}

ClusterSystem::ClusterAnalysis::ClusterAnalysis(const string path, const int iStep, const double accuracy)
{
	// read number of particles
	ifstream ifsRadius;
	ifsRadius.open(path + "output_radius.dat");
	if (ifsRadius.is_open()) {
		std::cout << "file is open\n";
	}
	int numParticle = 0;
	double dummyd;
	int dummyi;
	while (!ifsRadius.eof()) {
		ifsRadius >> dummyd;
		numParticle++;
	}
	numParticle--;
	std::cout << numParticle << std::endl;
	ifsRadius.close();

	// Generate particles from output files...
	_particles.reserve(numParticle);
	ifstream ifsRotation;
	ifsRotation.open(path + "output_rotation.dat");
	for (int i = 0; i < iStep - 1; ++i) {
		ifsRotation.ignore(2147483647, '\n');
	}
	for (int i = 0; i < numParticle; ++i) {
		ifsRotation >> dummyd;
		_particles.push_back(Particle(i, dummyd));
	}

	// Generate contacts from output file...
	ifstream ifsContact;
	ifsContact.open(path + "output_contact.dat");
	for (int i = 0; i < iStep - 1; ++i) {
		int numContact;
		ifsContact >> numContact;
		for (int j = 0; j < numContact; ++j) {
			/*ifsContact.ignore(1e10, '\n');*/
			int p1, p2;
			double fnx, fny, ftx, fty, fx, fy;
			ifsContact >> p1 >> p2 >> dummyd >> dummyd;
			ifsContact >> fnx >> fny >> ftx >> fty;
		}
	}
	int numContact;
	ifsContact >> numContact;
	_contacts.reserve(numContact);
	for (int i = 0; i < numContact; ++i) {
		int p1, p2;
		double fnx, fny, ftx, fty, fx, fy;
		double fn, ft, f;
		ifsContact >> p1 >> p2 >> dummyd >> dummyd;
		ifsContact >> fnx >> fny >> ftx >> fty;
		fx = fnx + ftx;
		fy = fny + fty;
		fn = sqrt(fnx * fnx + fny * fny);
		ft = sqrt(ftx * ftx + fty * fty);
		f = sqrt(fx * fx + fy * fy);
		_contacts.push_back(Contact(
			&_particles[p1],
			&_particles[p2],
			fn, ft, f));
		//std::cout << p1 << "\t" <<  p2 << std::endl;
		//std::cin.get();
	}

	// calculate threshold vector
	double dtheta = 0.0;
	while (dtheta < M_PI_2) {
		_dTheta0.push_back(dtheta);
		dtheta += accuracy;
	}
	_mclus.resize(_dTheta0.size(), vector<int>(_particles.size(), -1));
	_nclus.resize(_dTheta0.size(), vector<int>(_particles.size(), -1));
	_correlation.resize(_dTheta0.size(), 0.0);
}



void ClusterSystem::ClusterAnalysis::runAnalysis()
{
	calculateRelativeAngle();
#pragma omp parallel num_threads(8)
	{
#pragma omp for schedule(static)
		for (int t = 0; t < _dTheta0.size(); ++t) {
			vector<Cluster> clusters(_particles.size());
			for (int i = 0; i < clusters.size(); ++i) clusters[i] = Cluster(i);
			for (int i = 0; i < _contacts.size(); ++i) {
				Contact& c = _contacts[i];
				int ip = c._particle1->_id;
				int jp = c._particle2->_id;
				if (c._dTheta > _dTheta0[t]) continue;
				if (clusters[jp]._pnode == nullptr) {
					clusters[jp].linekTo(clusters[ip]);
				}
				else {
					if (clusters[ip].findHeadParentID() == clusters[jp].findHeadParentID()) continue;
					clusters[clusters[jp].findHeadParentID()].linekTo(clusters[ip]);
				}
			}
			std::vector<int> nclus(_particles.size(), -1);
			std::vector<int> mclus(_particles.size(), -1);
			for (int k1 = 0; k1 < _particles.size(); k1++) {
				nclus[k1] = k1;
				mclus[k1] = 1;
			}
			for (int i = 0; i < _particles.size(); ++i) {
				if (clusters[i]._pnode == nullptr) {

				}
				else {
					int index = clusters[i].findHeadParentID();
					nclus[i] = index;
					mclus[index]++;
					mclus[i] = 0;
				}
			}
			_mclus[t] = mclus;
			_nclus[t] = nclus;

			// remember this algo. it is the best way that I can do, to check if a contact is an internal contact inside a cluster
			// point biserial correlation coefficient: https://en.wikipedia.org/wiki/Point-biserial_correlation_coefficient
			vector<int> c_internal(_contacts.size(), 0);
			vector<int> c_boundary(_contacts.size(), 0);
			double count_internal = 0;
			double count_boundary = 0;
			double count_both = 0;
			double f_internal = 0.0;
			double f_boundary = 0.0;
			double f_both = 0.0;
			for (int i = 0; i < _contacts.size(); i++) {
				int gg1 = _contacts[i]._particle1->_id;
				int gg2 = _contacts[i]._particle2->_id;
				// if grain1 and grain2 are in the same cluster
				if (nclus[gg1] == nclus[gg2]) {
					c_internal[i] = 1;
					f_internal += _contacts[i]._totalForce;
					count_internal++;
				}
				else {
				//else if (mclus[nclus[gg1]] > 1 || mclus[nclus[gg2]] > 1) {
					c_boundary[i] = 1;
					f_boundary += _contacts[i]._totalForce;
					count_boundary++;
				}
			}
			f_both = f_internal + f_boundary;
			count_both = count_internal + count_boundary; // n = n1 + n0
			f_internal /= count_internal; // M1
			f_boundary /= count_boundary; // M0
			f_both /= count_both; // x_bar
			double sn = 0.0; // sn
			for (int i = 0; i < _contacts.size(); i++) {
				if (c_internal[i] == 1 || c_boundary[i] == 1) {
					double temp = _contacts[i]._totalForce - f_both;
					sn += temp * temp;
				}
			}
			sn /= count_both;
			sn = std::sqrt(sn);
			double rpb = (f_internal - f_boundary) / sn * std::sqrt(count_internal * count_boundary / count_both / count_both);
			for (int i = 0; i < _contacts.size(); i++) {
				if (c_internal[i] == 1 && c_boundary[i] == 1) std::cout << "Error: c_internal = c_boundary, i = " << i << std::endl;
			}

			if (isnormal(rpb)) _correlation[t] = rpb;
			
		}
	}
}

void ClusterSystem::ClusterAnalysis::printClusters2File(string path)
{
	
	/*stringstream md;
	md << "md " << path.c_str();
	std::cout << md.str() << std::endl;
	system(md.str().c_str());*/

	std::stringstream fileName;
	fileName << path << "clusterLabelling_mclus.dat";
	FILE* fileMclus = fopen(fileName.str().c_str(), "w");
	fileName.str("");
	fileName << path << "clusterLabelling_nclus.dat";
	FILE* fileNclus = fopen(fileName.str().c_str(), "w");
	fileName.str("");
	fileName << path << "clusterLabelling_correlation.dat";
	FILE* fileCoorelation = fopen(fileName.str().c_str(), "w");

	for (int i = 0; i < _dTheta0.size(); ++i) {
		for (int j = 0; j < _particles.size(); ++j) {
			fprintf(fileMclus, "%d\t", _mclus[i][j]);
			fprintf(fileNclus, "%d\t", _nclus[i][j]);
		}
		fprintf(fileMclus, "\n");
		fprintf(fileNclus, "\n");
		fprintf(fileCoorelation, "%le\t%le\n", _dTheta0[i], _correlation[i]);
	}
	fclose(fileMclus);
	fclose(fileNclus);
	fclose(fileCoorelation);
}

void ClusterSystem::ClusterAnalysis::calculateRelativeAngle()
{
#pragma omp parallel num_threads(4)
	{
#pragma omp for schedule(static)
		for (int i = 0; i < _contacts.size(); ++i) {
			_contacts[i].calculateRelativeAngle();
		}
	}
}

ClusterSystem::Cluster::Cluster(int id) : _id(id)
{

}

int ClusterSystem::Cluster::findHeadParentID()
{
	Cluster* headParent = this;
	while (headParent->_pnode != nullptr) {
		headParent = headParent->_pnode;
	}
	return headParent->_id;
}

void ClusterSystem::Cluster::linekTo(Cluster& parent)
{
	if (_pnode != nullptr) {
		std::cout << "__ERROR__this cluster already has a parent cluster...\n";
		exit(1);
	}
	_pnode = &parent;
	return;
}
