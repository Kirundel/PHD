using namespace std;

#define _CRT_SECURE_NO_WARNINGS
#include <vector>
#include <map>
#include <algorithm>
#include <iomanip>
#include <random>
#include <chrono>

//#include <fstream>
//ifstream cin("input.txt");
//ofstream cout("output.txt");
#include <iostream>


using ll = long long;
using pii = pair<int, int>;
using pdd = pair<double, double>;

template <typename T>
inline void sort(vector<T>& x)
{
	sort(x.begin(), x.end());
}

template <typename T>
inline void unique(vector<T>& x)
{
	sort(x);
	auto end_s = unique(x.begin(), x.end()) - x.begin();
	x.resize(end_s);
}

struct fenwick_tree
{
	vector<ll> a;

	fenwick_tree(int n) : a(n + 2)
	{}

	void add(int pos, ll delta)
	{
		for (pos++; pos < a.size(); pos += pos & -pos)
			a[pos] += delta;
	}

	ll get_sum(int pos)
	{
		ll sum = 0;
		for (pos++; pos > 0; pos -= pos & -pos)
			sum += a[pos];
		return sum;
	}

	ll get_sum(int l, int r)
	{
		return this->get_sum(r) - this->get_sum(l - 1);
	}
};

class compression
{
public:
	vector<double> a;

	void clear() {
		a.clear();
	}

	void add_value(double v)
	{
		a.push_back(v);
	}

	void compress()
	{
		unique(a);
	}

	int get_val(double value)
	{
		return lower_bound(a.begin(), a.end(), value) - a.begin();
	}

	int get_up_val(double value)
	{
		return upper_bound(a.begin(), a.end(), value) - a.begin();
	}
};

int n, m;
ll k;
double precision;
vector<pdd> points;
compression b_compression;

using iteration_pairs = pair<double, pair<int, double>>;

map<int, vector<iteration_pairs>>vectors_merge_map;

vector<iteration_pairs>& merge(vector<iteration_pairs>& first_vector, vector<iteration_pairs>& second_vector)
{
	vector<iteration_pairs>& result_vector = vectors_merge_map[first_vector.size() + second_vector.size()];

	int l = 0, r = 0;
	while (l < first_vector.size() || r < second_vector.size())
	{
		if (l >= first_vector.size())
		{
			result_vector[l + r] = second_vector[r];
			r++;
		}
		else
		{
			if (r >= second_vector.size())
			{
				result_vector[l + r] = first_vector[l];
				l++;
			}
			else
			{
				if (first_vector[l] < second_vector[r])
				{
					result_vector[l + r] = first_vector[l];
					l++;
				}
				else
				{
					result_vector[l + r] = second_vector[r];
					r++;
				}
			}
		}
	}
	return result_vector;
}

vector<iteration_pairs> events_0;
vector<iteration_pairs> events_1;
vector<iteration_pairs> events_2;


ll count_with_length(double deltas)
{
	for (int i = 0; i < points.size(); i++)
	{
		auto& p = points[i];
		events_0[i] = { p.first - deltas, { 0, p.second } };
		events_1[i] = { p.first + deltas + 1, { 1 , p.second } };
		events_2[i] = { p.first, { 2 , p.second } };
	}

	auto& events_t = merge(events_0, events_1);
	auto& events = merge(events_t, events_2);

	fenwick_tree fwt(b_compression.a.size() + 3);

	ll answer = 0;

	for (auto ev : events)
	{
		auto event_y_coord = ev.second.second;
		int event_type = ev.second.first;
		if (event_type == 2)
		{
			int left = b_compression.get_val(event_y_coord - deltas);
			int right = b_compression.get_up_val(event_y_coord + deltas) - 1;
			answer += fwt.get_sum(left, right) - 1;
		}
		else
		{
			event_y_coord = b_compression.get_val(event_y_coord);
			if (event_type == 0)
			{
				fwt.add(event_y_coord, 1);
			}
			if (event_type == 1)
			{
				fwt.add(event_y_coord, -1);
			}
		}

	}

	return answer / 2;
}


vector<pair<double, double>> _data;

void init_data(int _n, int _k, double _precision) {
	n = _n;
	k = _k;
	precision = _precision;
	_data.resize(_n);

	std::uniform_real_distribution<double> unif(0.0, 1.0);
	std::default_random_engine re;
	
	for (int i = 0; i < n; i++) {
		_data[i] = { unif(re), unif(re) };
	}

}

void init()
{
	vectors_merge_map[2 * n] = vector<iteration_pairs>(2 * n);
	vectors_merge_map[3 * n] = vector<iteration_pairs>(3 * n);
	events_0.assign(n, { 0.0, {0, 0.0} });
	events_1.assign(n, { 0.0, {0, 0.0} });
	events_2.assign(n, { 0.0, {0, 0.0} });

}

double solve()
{
	init();

	for (int i = 0; i < n; i++)
	{
		double x = _data[i].first, y = _data[i].second;
		double a = x + y;
		double b = y - x;
		points.push_back({ a, b });
		b_compression.add_value(b);
	}

	sort(points);

	b_compression.compress();

	double left = 0.0;
	double right = 4;

	while (right - left > precision)
	{
		double mid = (left + right) / 2;
		if (count_with_length(mid) < k)
		{
			left = mid;
		}
		else
		{
			right = mid;
		}
	}

	auto result = (left + right) / 2;
	return result;
}

void clear() {
	points.clear();
	b_compression.clear();
	points.reserve(n);
	vectors_merge_map.clear();
}

void calculate_time(int _n, int _k, double _precision) {
	init_data(_n, _k, _precision);

	auto result = 0.0;
	ll sum_time = 0;

	cout << _n << " " << _k << " " <<  " : ";
	int number_runs = 1;
	for (int i = 0; i < number_runs; i++) {
		clear();
		auto begin = std::chrono::steady_clock::now();
		result += solve();
		auto end = std::chrono::steady_clock::now();
		sum_time += (std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count());
	}

	if (result > 0.000000001) {
		cout << (ll)(sum_time / (1e6 * number_runs) + 0.5) << endl;
	}
}

int main()
{
	for (int i = 16; i <= 25; i++) {
		calculate_time(1 << i, 1000, 1e-5);
	}

	for (int i = 17; i <= 25; i++) {
		calculate_time(1 << i, 100000, 1e-5);
	}

	return 0;
}