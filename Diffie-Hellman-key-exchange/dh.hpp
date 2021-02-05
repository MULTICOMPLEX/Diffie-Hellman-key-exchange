
#include <cstdio>
#include <chrono>
#include <iostream>
#include <iomanip>
#include <unordered_map>
#include <algorithm>
#include <vector>
#include <tuple>
#include <set>

#include <boost/multiprecision/cpp_int.hpp>
#include <boost/multiprecision/miller_rabin.hpp>
#include <boost/integer/mod_inverse.hpp>
#include <boost/lexical_cast.hpp>
using namespace boost::multiprecision;
using namespace boost::multiprecision::literals;

typedef int1024_t INT;

#include "mxws.hpp"

template<typename T>
class mxws_t
{

private:

	std::random_device r;

public:

	static const unsigned n = std::numeric_limits<T>::digits * 2;
	typedef number<cpp_int_backend<n * 2, n * 2, unsigned_magnitude, unchecked, void> > uint;

	uint x, w;
	
	typedef T result_type;

	void seed()
	{
		init();
		x = 1;
	}

	mxws_t(const std::seed_seq& seq)
	{	
		if (seq.size() * 32 == n)
		{
			std::vector<uint32_t> seeds(seq.size());
			seq.param(seeds.begin());
			for (unsigned x = 0; x < n; x += 32)
				w |= (uint(seeds[x / 32]) << x);
			x = 1;
		}

		else init();
	}

	mxws_t(const uint& seed)
	{
		for (unsigned x = 0; x < n; x += 32)
			w |= (uint(r()) << x);
		x = 1;
	}

	void init(const uint64_t& seed)
	{
		for (unsigned x = 0; x < n; x += 32)
			w |= (uint(r()) << x);
		x = 1;
	}

	mxws_t()
	{
		init();
	}

	void init()
	{
		for (unsigned x = 0; x < n; x += 32)
			w |= (uint(r()) << x);
		x = 1;
	}

	virtual ~mxws_t() = default;

	T min() { return std::numeric_limits<T>::min(); }
	T max() { return std::numeric_limits<T>::max(); }

	inline T operator()()
	{
		x *= w;
		x = (x >> (n/2)) | (x << (n/2));
		w += x;

		return T(x);
	}

	inline T operator()(const T& max)
	{
		return (*this)() % (max + 1);
	}

	inline T operator()(const T& min, const T& max)
	{
		return min + ((*this)() % (max - min + 1));
	}

};

//////////////////////////////////////////

inline uint64_t rand_uint64()
{
	static thread_local mxws_64 mxws_64;
	return mxws_64();
}

//////////////////////////////////////////

template<typename T>
inline T rand_t_prime()
{

	static thread_local mxws_t<T> mxws_t;
	T n;
	while (!miller_rabin_test(n, 1, mxws_t))n = mxws_t();
	return n;
}

template<typename T>
inline T rand_t()
{
	static thread_local mxws_t<T> mxws_t;
	return mxws_t();
}

//////////////////////////////////////////

template<typename T>
inline bool isPrime
(
	const T& n
)
{

	static thread_local mxws_t<T> mxws_t;

	if (!miller_rabin_test(n, 10, mxws_t))
		return false;

	return true;
}

template<typename T>
inline T Primitive_Root
(
	const T& p
)
{
	T r;
	
	mxws_t<T> mxws_t;

	while (!((gcd(r, p) == 1) && isPrime<T>(r)))r = mxws_t();

	return r;
}

const thread_local uint512_t P = rand_t_prime<uint512_t>();
const thread_local uint512_t G = Primitive_Root(P);

template<typename T>
inline void dh_gen_keypair
(
	T& private_key,
	T& public_key
)
{
	private_key = rand_t<T>();
	public_key = powm(G, private_key, P);
}

template<typename T>
inline T dh_gen_secret
(
	T my_private_key,
	T another_public_key
)
{
	return powm(another_public_key, my_private_key, P);
}

/*
		An element in a Galois field FP
		Adapted for the specific behaviour of the "mod" function where (-n) mod m returns a negative number
		Allows basic arithmetic operations between elements:
				+,-,/,scalar multiply

		The template argument P is the order of the field
*/

template<typename T>
class FiniteFieldElement
{

public:

	T i_;

	static T Prime;

	FiniteFieldElement() = default;
	virtual ~FiniteFieldElement() = default;

	inline void assign(const T& i)
	{
		i_ = i;
		if (i < 0)
		{
			// ensure (-i) mod p correct behaviour
			// the (i%P) term is to ensure that i is in the correct range before normalizing
			i_ = (i % Prime) + 2 * Prime;
		}

		i_ %= Prime;
	}

	// ctor
	explicit FiniteFieldElement(T i)
	{
		assign(i);
	}

	// copy ctor
	FiniteFieldElement(const FiniteFieldElement& rhs)
		: i_(rhs.i_) {}

	// access "raw" integer
	T i() const { return i_; }

	// negate
	inline FiniteFieldElement operator-() const
	{
		return FiniteFieldElement(-i_);
	}

	// assign from integer
	inline FiniteFieldElement operator=(T i)
	{
		assign(i);
		return *this;
	}

	// assign from field element
	inline FiniteFieldElement operator=(const FiniteFieldElement& rhs)
	{
		i_ = rhs.i_;
		return *this;
	}

	// *=
	inline FiniteFieldElement operator*=(const FiniteFieldElement& rhs)
	{
		i_ = (i_ * rhs.i_) % P;
		return *this;
	}

	// ==
	inline friend bool operator==(const FiniteFieldElement& lhs, const FiniteFieldElement& rhs)
	{
		return (lhs.i_ == rhs.i_);
	}

	// == int
	inline friend bool operator==(const FiniteFieldElement& lhs, T rhs)
	{
		return (lhs.i_ == rhs);
	}

	// !=
	inline friend bool operator!=(const FiniteFieldElement& lhs, T rhs)
	{
		return (lhs.i_ != rhs);
	}

	// a / b
	friend FiniteFieldElement operator/(const FiniteFieldElement& lhs, const FiniteFieldElement& rhs)
	{
		return FiniteFieldElement(lhs.i_ * boost::integer::mod_inverse(rhs.i_, Prime));
	}

	// a + b
	inline friend FiniteFieldElement operator+(const FiniteFieldElement& lhs, const FiniteFieldElement& rhs)
	{
		return FiniteFieldElement(lhs.i_ + rhs.i_);
	}

	// a - b
	inline friend FiniteFieldElement operator-(const FiniteFieldElement& lhs, const FiniteFieldElement& rhs)
	{
		return FiniteFieldElement(lhs.i_ - rhs.i_);
	}

	// a + int
	inline friend FiniteFieldElement operator+(const FiniteFieldElement& lhs, T i)
	{
		return FiniteFieldElement(lhs.i_ + i);
	}

	// int + a
	inline friend FiniteFieldElement operator+(T i, const FiniteFieldElement& rhs)
	{
		return FiniteFieldElement(rhs.i_ + i);
	}

	// int * a
	inline friend FiniteFieldElement operator*(T n, const FiniteFieldElement& rhs)
	{
		return FiniteFieldElement(n * rhs.i_);
	}

	// a * b
	friend FiniteFieldElement operator*(const FiniteFieldElement& lhs, const FiniteFieldElement& rhs)
	{
		return FiniteFieldElement(lhs.i_ * rhs.i_);
	}

	// ostream handler
	inline friend  std::ostream& operator<<(std::ostream& os, const FiniteFieldElement& g)
	{
		return os << g.i_;
	}
};


/*
		Elliptic Curve over a finite field of order P:
		y^2 mod P = x^3 + ax + b mod P

		Template parameter P is the order of the finite field Fp over which this curve is defined
*/

template<typename T>
class EllipticCurve : public FiniteFieldElement<T>
{

public:
	// this curve is defined over the finite field (Galois field) Fp, this is the 
	// typedef of elements in it
	typedef FiniteFieldElement<T> ffe_t;

	EllipticCurve() = default;
	virtual ~EllipticCurve() = default;

	/*
			A point, or group element, on the EC, consisting of two elements of the field FP
			Points can only created by the EC instance itself as they have to be
			elements of the group generated by the EC
	*/

	class Point
	{
	public:

		ffe_t  x_;
		ffe_t  y_;
		EllipticCurve* ec_;

		Point() = default;
		virtual ~Point() = default;

		// core of the doubling multiplier algorithm (see below)
		// multiplies acc by m as a series of "2*acc's"

		void addDouble(const T& m, Point& acc)
		{
			if (m > 0)
			{
				Point r = acc;
				for (T n = 0; n < m; ++n)
				{
					r += r;     // doubling step                          
				}
				acc = r;
			}
		}

		// doubling multiplier algorithm
		// multiplies a by k by expanding in multiplies by 2
		// a is also an accumulator that stores the intermediate results
		// between the "1s" of the binary form of the input scalar k
		Point scalarMultiply(const T& k, const Point& a)
		{
			Point acc = a;
			Point res = Point(0, 0, *ec_);
			T i = 0, j = 0;
			T b = k;

			while (b)
			{
				if (b & 1)
				{
					// bit is set; acc = 2^(i-j)*acc
					addDouble(i - j, acc);
					res += acc;
					j = i;  // last bit set
				}
				b >>= 1;
				++i;
			}
			return res;
		}

		// adding two points on the curve
		void add(ffe_t x1, ffe_t y1, ffe_t x2, ffe_t y2, ffe_t& xR, ffe_t& yR)
		{
			// special cases involving the additive identity                     
			if (x1 == 0 && y1 == 0)
			{
				xR = x2;
				yR = y2;
				return;
			}
			if (x2 == 0 && y2 == 0)
			{
				xR = x1;
				yR = y1;
				return;
			}
			if (y1 == -y2)
			{
				xR = yR = 0;
				return;
			}

			// the additions
			ffe_t s;

			if (x1 == x2 && y1 == y2)
			{
				//2P                          
				s = (3 * (x1.i() * x1.i()) + ec_->a()) / (2 * y1);
				xR = ((s * s) - 2 * x1);
			}
			else
			{
				//P+Q                            
				s = (y1 - y2) / (x1 - x2);
				xR = ((s * s) - x1 - x2);
			}

			if (s != 0)
			{
				yR = (-y1 + s * (x1 - xR));
			}
			else
			{
				xR = yR = 0;
			}
		}

		Point(const T& x, const T& y)
			: x_(x),
			y_(y),
			ec_(0)
		{}

		Point(const T& x, const T& y, EllipticCurve& EllipticCurve)
			: x_(x),
			y_(y),
			ec_(&EllipticCurve)
		{}

		Point(const ffe_t& x, const ffe_t& y, EllipticCurve& EllipticCurve)
			: x_(x),
			y_(y),
			ec_(&EllipticCurve)
		{}

		// copy ctor
		Point(const Point& rhs)
		{
			x_ = rhs.x_;
			y_ = rhs.y_;
			ec_ = rhs.ec_;
		}

		// assignment
		Point operator=(const Point& rhs)
		{
			x_ = rhs.x_;
			y_ = rhs.y_;
			ec_ = rhs.ec_;
			return *this;
		}

		// access x component as element of Fp
		ffe_t x() const { return x_; }
		// access y component as element of Fp
		ffe_t y() const { return y_; }

		// negate
		Point operator-()
		{
			return Point(x_, -y_);
		}

		// ==
		bool operator== (const Point& rhs)
		{
			return (*this.ec_ == rhs.ec_) && (*this.x_ == rhs.x_) && (*this.y_ == rhs.y_);
		}

		// !=
		bool operator!=(const Point& rhs)
		{
			return (*this.ec_ != rhs.ec_) || (*this.x_ != rhs.x_) || (*this.y_ != rhs.y_);
		}

		// a + b         
		Point operator+(const Point& rhs) const
		{
			ffe_t xR, yR;
			*this.add(*this.x_, *this.y_, rhs.x_, rhs.y_, xR, yR);
			return *this(xR, yR, *this.ec_);
		}

		// a * int
		friend Point operator*(const T& k, const Point& rhs)
		{
			return Point(rhs).operator*=(k);
		}

		// +=
		Point operator+=(const Point& rhs)
		{
			add(x_, y_, rhs.x_, rhs.y_, x_, y_);
			return *this;
		}

		// a *= int
		Point operator*=(const T& k)
		{
			return (*this = scalarMultiply(k, *this));
		}

		// ostream handler: print this point
		friend std::ostream& operator <<(std::ostream& os, const Point& p)
		{
			return (os << "(" << p.x_ << ", " << p.y_ << ")");
		}
	};

	// ==================================================== EllipticCurve impl

	// ctor
	// Initialize EC as y^2 = x^3 + ax + b
	EllipticCurve(const T& a, const T& b)
		: a_(a),
		b_(b)
	{}

	const Point G1 = Point(
		T("0xa169f51dbe7886b50979ee4568b1d2a5432e23c2459f159b46bed404b79e38b8"),
		T("0xc28263f9c1f0d94478f908406817638b31786d464ec1bca932c90471b5ead9d3"), *this);

	// the parameter a (as an element of Fp)
	FiniteFieldElement<T> a() const { return a_; }
	// the paramter b (as an element of Fp)
	FiniteFieldElement<T> b() const { return b_; }

	// ostream handler: print this curve in human readable form
	template <typename A>
	friend std::ostream& operator <<(std::ostream& os, const EllipticCurve<A>& EllipticCurve);

	FiniteFieldElement<T>       a_;         // paramter a of the EC equation
	FiniteFieldElement<T>       b_;         // parameter b of the EC equation

};

template <typename A>
std::ostream& operator <<(std::ostream& os, const EllipticCurve<A>& EllipticCurve)
{
	os << "y^2 mod " << EllipticCurve.Prime << " = (x^3 ";

	std::string s1, s2;
	if (EllipticCurve.a_.i_ < 0)
		s1 = "-";
	else s1 = "+";

	if (EllipticCurve.b_.i_ < 0)
		s2 = "-";
	else s2 = "+";

	if (EllipticCurve.a_ != 0)
	{
		os << s1 << EllipticCurve.a_ << "x ";
	}

	if (EllipticCurve.b_.i() != 0)
	{
		os << s2 << EllipticCurve.b_;
	}

	os << ") mod " << std::endl << EllipticCurve.Prime;
	return os;
}

template <typename T>
class Solution
{
public:
	T root1, root2;
	bool exists;
};

template <typename T>
class Solution <T> makeSolution(const T& root1, const T& root2, bool exists)
{
	class Solution <T> sol;
	sol.root1 = root1;
	sol.root2 = root2;
	sol.exists = exists;
	return sol;
}

template <typename T>
class Solution<T> ts(const T& n, const T& p)
{
	T q = p - 1;
	T ss = 0;
	T z = 2;
	T c, r, t, m;

	if (powm(n, (p - 1) / 2, p) != 1) {
		return makeSolution(T(0), T(0), false);
	}

	while ((q & 1) == 0) {
		ss += 1;
		q >>= 1;
	}

	if (ss == 1) {
		T r1 = powm(n, (p + 1) / 4, p);
		return makeSolution<T>(r1, p - r1, true);
	}

	while (powm(z, (p - 1) / 2, p) != p - 1) {
		z++;
	}

	c = powm(z, q, p);
	r = powm(n, (q + 1) / 2, p);
	t = powm(n, q, p);
	m = ss;

	while (true) {
		T i = 0, zz = t;
		T b = c, e;

		if (t == 1) {
			return makeSolution<T>(r, p - r, true);
		}

		while (zz != 1 && i < (m - 1)) {
			zz = zz * zz % p;
			i++;
		}
		e = m - i - 1;

		while (e > 0) {
			b = b * b % p;
			e--;
		}

		r = r * b % p;
		c = b * b % p;
		t = t * c % p;
		m = i;
	}
}

template<typename T>
T compress(const T& PubKey)
{
	int e;
	if ((PubKey.y().i_ & 1)) e = 1;
	else e = 0;

	const T v(PubKey.x().i_, INT(e));
	return v;
}

template<typename T, typename P, typename Pr>
T decompress(const T& PubKey, const P& Prime, const Pr& EllipticCurve)
{
	//sqrtmod(pow(x, 3, p) + a * x + b, p)

	class Solution<P> sol =
		ts(powm(PubKey.x().i_, 3, Prime) + (EllipticCurve.a_.i() * PubKey.x().i_) + EllipticCurve.b_.i(), Prime);

	const T pointa(PubKey.x().i_, sol.root1);
	const T pointb(PubKey.x().i_, sol.root2);

	auto oe = int(PubKey.y().i_ & 1);

	auto oe_r1 = int(sol.root1 & 1);

	if (oe && oe_r1)
	{
		return pointa;
	}

	if (!oe && !oe_r1)
	{
		return pointa;
	}

	return pointb;
}

