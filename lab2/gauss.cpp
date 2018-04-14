#include "gauss.h"

#include <cmath>

const double inf = 100000000.0, eps = 0.00000001;

int gauss(vector <double> &x, vector <vector <double>> &a, vector <double> &b) {
	x.clear();

	if (a.size() != b.size()) {
		return 1;
	}
	for (const auto &i : a) {
		if (i.size() != b.size()) {
			return 1;
		}
	}

	vector <vector <double>> tmpa;
	vector <double> tmpb;
	vector <int> tmpflag;

	for (const auto &i : a) {
		tmpa.emplace_back(0);
		for (const auto &j : i) {
			tmpa.back().emplace_back(j);
		}
	}

	for (const auto &i : b) {
		tmpb.emplace_back(i);
	}

	double tmpmax, tmpr;
	int tmpp;

	for (int i = 0; i < a.size(); i++) {
		tmpp = -1;
		tmpmax = -inf;
		for (int j = 0; j < a[i].size(); j++) {
			if (!tmpflag[j] && abs(tmpa[j][i]) > tmpmax) {
				tmpp = j;
				tmpmax = abs(tmpa[j][i]);
			}
		}

		if (tmpp == -1) {
			return 1;
		}
		
		tmpr = tmpa[tmpp][i];
		for (int j = 0; j < tmpa[i].size(); j++) {
			tmpa[tmpp][j] /= tmpr;
		}
		tmpb[tmpp] /= tmpr;
		tmpflag[tmpp] = 1;

		for (int j = 0; j < tmpa.size(); j++) {
			if (tmpp != j) {
				tmpr = tmpa[j][i];
				for (int k = 0; k < tmpa[i].size(); k++) {
					tmpa[j][k] -= tmpr * tmpa[tmpp][k];
				}
				tmpb[j] -= tmpr * tmpb[tmpp];
			}
		}
	}
	for (int i = 0; i < tmpa.size(); i++) {
		x.emplace_back(0);
	}
	for (int i = 0; i < tmpa.size(); i++) {
		for (int j = 0; j < tmpa[i].size(); j++) {
			if (tmpa[i][j] > eps) {
				x[j] = tmpb[j];
			}
		}
	}
	return 0;
}
