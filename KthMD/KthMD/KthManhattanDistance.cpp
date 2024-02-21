using namespace std;

#define _CRT_SECURE_NO_WARNINGS
#include <vector>
#include <map>
#include <algorithm>

#include <fstream>
ifstream cin("input.txt");
ofstream cout("output.txt");


using ll = long long;
using pii = pair<int, int>;

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
	vector<int> a;

	void add_value(ll v)
	{
		a.push_back(v);
	}

	void compress()
	{
		unique(a);
	}

	int get_val(int value)
	{
		return lower_bound(a.begin(), a.end(), value) - a.begin();
	}

	int get_up_val(int value)
	{
		return upper_bound(a.begin(), a.end(), value) - a.begin();
	}
};

int n, m;
ll k;
vector<pii> points;
compression b_compression;

map<int, vector<pair<int, pii>>>vectors_merge_map;

template<typename T>
vector<T>& merge(vector<T>& first_vector, vector<T>& second_vector)
{
	vector<T>& result_vector = vectors_merge_map[first_vector.size() + second_vector.size()];

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

vector<pair<int, pii>> events_0;
vector<pair<int, pii>> events_1;
vector<pair<int, pii>> events_2;


ll count_with_length(int deltas)
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
		int event_y_coord = ev.second.second;
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

void init()
{
	vectors_merge_map[2 * n] = vector<pair<int, pii>>(2 * n);
	vectors_merge_map[3 * n] = vector<pair<int, pii>>(3 * n);
	events_0.resize(n);
	events_1.resize(n);
	events_2.resize(n);
}

void solve()
{
	cin >> n >> k;

	init();

	for (int i = 0; i < n; i++)
	{
		int x, y;
		cin >> x >> y;
		int a = x + y;
		int b = y - x;
		points.push_back({ a, b });
		b_compression.add_value(b);
	}

	sort(points);

	b_compression.compress();

	int left = -1;
	int right = 4e8 + 1;

	while (right - left > 1)
	{
		int mid = (left + right) / 2;
		if (count_with_length(mid) < k)
		{
			left = mid;
		}
		else
		{
			right = mid;
		}
	}

	cout << right << endl;
}

int main()
{
	solve();

	return 0;
}