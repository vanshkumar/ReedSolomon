import random
import itertools

# Field of size p = q^m, (n = q-1, k) code that is t error-correcting
p = 47
m = 1
q = pow(p, m)
n = q - 1
t = 4
k = n - 2 * t

# Primitive element that generates finite field of size p
alpha = 5

# Create the exponential table by computing powers of alpha
exp = [1]
last = alpha
while exp[-1] != 1 or len(exp) < 2:
	exp = exp + [last]
	last = (last * alpha) % q

# Create the log table from the exponential table 
logs = [None]
for i in range(1, q):
	logs.append(exp.index(i))

# A finite field integer class that abstracts away all finite field arithmetic
# Note that it subclasses int, as every element of a prime field can be mapped
# to an integer
class GFInt(int):
	# Initialize the only attribute of this class, an integer representation
	def __init__(self, value):
		self.integer = value % q

	# Addition modulo q (field size)
	def __add__(self, other):
		return GFInt((self.integer + other.integer) % q)

	# Extension of addition
	def __sub__(self, other):
		return self + (GFInt(-1) * other)

	# Multiplying by using exponential and log lookup tables. Multiplication
	# just becomes exponentiation of the sum of the logs of the two integers
	# we are multiplpying, mod q-1 of course. It is q-1 instead of q because
	# this is in the multiplicative group
	def __mul__(self, other):
		myPower = logs[self.integer]
		otherPower = logs[other.integer]

		# If either is 0, return 0
		if myPower == None or otherPower == None:
			return GFInt(0)
		return GFInt(exp[(myPower + otherPower) % (q - 1)])

	# Similar to multiplication in that use lookup tables; subtract logs of the
	# two integers instead of adding them.
	def __div__(self, other):
		myPower = logs[self.integer]
		otherPower = logs[other.integer]

		# If the numerator is 0, return 0
		if myPower == None:
			return GFInt(0)

		# If the denominator is 0, return None to signify a divide by 0 error
		if otherPower == None:
			#print "Divide by zero error in div!"
			return None
		return GFInt(exp[(myPower - otherPower) % (q - 1)])


	# Comparison testing, very similar to just that for integers 

	def __eq__(self, other):
		if isinstance(other, GFInt):
			return self.integer == other.integer
		elif isinstance(other, int):
			return self.integer == other

	def __ne__(self, other):
		return not (self == other)

	def __gt__(self, other):
		if isinstance(other, GFInt):
			return self.integer > other.integer
		elif isinstance(other, int):
			return self.integer > other

	def __lt__(self, other):
		return not (self > other) and (self != other)

	# Exponentiate by returning the exponential of the log of the integer
	# times m (what we want to take self to the power of) mod q-1
	def exponentiate(self, m):
		return GFInt(exp[(logs[self.integer] * m) % (q - 1)])

	# Inverse is a lookup of the log table, just the exponential of
	# q - 1 - the log of the integer
	def inverse(self):
		return GFInt(exp[q - 1 - logs[self.integer]])

# Add two polynomials, a and b. This just boils down to adding element-wise, as
# the (mod q) is abstracted away as all elements are GFInt's
def polyAdd(a, b):
	# Make the two the same length
	if len(a) < len(b):
		a = a + [GFInt(0)]*(len(b) - len(a))
	elif len(b) < len(a):
		b = b + [GFInt(0)]*(len(a) - len(b))

	return [a[i] + b[i] for i in range(len(a))]

# Perform a - b, extension of polynomial addition
def polySub(a, b):
	return polyAdd(a, [GFInt(-1) * x for x in b])

# Multiply two polynomials by just looping through the terms of one of them,
# computing the product of that term with the other polynomial, and summing
# all of these computations
def polyMult(a, b):
	prod = [GFInt(0)]*(len(a) + len(b) - 1)
	for i in range(len(a)):
		tmp = [GFInt(0)] * i + [x * a[i] for x in b]
		prod = polyAdd(prod, tmp)
	return prod

# Return the degree of a given polynomial. Boils down to finding the first
# nonzero coefficient of the polynomial, starting from the end (highest power)
def degree(poly):
	# Convert to binary, where 1 means nonzero coefficient
	check = [int(x > 0) for x in poly]

	# If the polynomial is 0, return -1 to signify this
	if check.count(GFInt(1)) == 0:
		return -1

	# Reverse because we want the first instance of a 1
	check.reverse()
	return len(check) - 1 - check.index(GFInt(1))

# Divide a by b. Uses long divisions
def polyDiv(a, b):
	# If a is 0, return [0, 0]
	if degree(a) < 0:
		return [GFPoly([0] * len(a)), GFPoly([0] * len(a))]

	# If b is 0, we are trying to divide by zero
	if degree(b) < 0:
		print "Divide by zero error"

	# Initialize quotient and remainder
	quot = [GFInt(0)] * len(a)
	rem  = a

	# If degree(b) > degree(a), then the quotient is just 0 and the remainder
	# is simply b
	if degree(a) >= degree(b):

		# While we can still divide, do so
		x = degree(a) - degree(b)
		while x >= 0:
			# Divisor at a given stage of the calculation
			div = [GFInt(0)] * x  + b

			# Set next value of the quotient
			quot[x] = a[degree(a)] / div[degree(div)]

			# Calculate what we need to subtract from a and subtract from it
			prod = [d * quot[x] for d in div]
			a = polySub(a, prod)

			# New difference in degree of the two
			x = degree(a) - degree(b)

		# Remainder is whatever's left over
		rem = a
	return [quot, rem]

# Convert a list of integers to a list of GFInts so operations can be performed
def GFPoly(lst):
	return [GFInt(x) for x in lst]

# Compute the generator polynomial. This just multiplies out
# (1-alpha)(1-alpha^2)...(1-alpha^(2t)
def computeGenerator():
	gener = [GFInt(1)]
	for i in range(1, n-k+1):
		tmp = [GFInt(-1*exp[i]), GFInt(1)]
		gener = polyMult(gener, tmp)
	return gener

# Encode a message m(x) as c(x) = x^(2t)*m(x) - r(x), where r(x) is the
# remainder when x^(2t)*m(x) is divided by g(x), the generator
def encode(m):
	# Multiply b x^(2t) by adding 2t zeros to the left of m
	tmp = [0] * (2*t) + m

	# Pad with zeros on the right if necessary
	shift = tmp + [GFInt(0)] * int(len(gener) >= len(tmp)) * (len(gener) - len(tmp))
	shift = GFPoly(shift)
	shifted = shift[:]

	# Divide and get remainder
	[quot, rem] = polyDiv(shift, gener)

	return [shifted[i] - rem[i] for i in range(len(rem))]

# Send a codeword through the channel with probability prob
def channel(c, prob):
	# Error vector
	e = [GFInt(0)] * len(c)

	# Table of probabilities, last value of this should be 1
	probs = [GFInt(0)] * (q+1)

	# Probability of having error value be 0
	probs[0] = 1 - prob

	# Probability of having it be nonzero
	nonzero = prob / (q - 1)

	# Generate table/bins
	for i in range(1, q+1):
		probs[i] = probs[i-1] + nonzero

	# For each error location, compute random number and see which bin it falls
	# into in the table of probabilities. That is e_i.
	for i in range(len(e)):
		p = random.random()
		check = [int(p < x) for x in probs]
		e[i] = GFInt(check.index(1))

	# Return y = c + e
	return polyAdd(c, e)

# Compute the syndrome of a received vector y by evaluating at various powers
# of alpha, the primitive element
def syndrome(y):
	syn = [GFInt(0)] * (2 * t)
	for i in range(1, 2 * t + 1):
		tmp = [y[j] * GFInt(alpha).exponentiate(i * j) for j in range(len(y))]
		syn[i-1] = GFInt(sum(tmp) % q)
	return syn

# Compute the error locator and error evaluator polynomials via the extended
# Euclidean algorithm
def computeErrors(syn):
	# Intialize: 
		# R is remainder, S is syndrome, A is coefficient polynomial of S, B is
		# coefficient polynomial of R, and i is counter of iterations
	R = GFPoly([0] * (2*t) + [1])
	S0 = syn
	S = S0[:]
	A = [GFInt(1)]
	B = [GFInt(0)]
	i = 0

	# Loop until we cannot divide anymore
	while degree(S) >= t:
		# Temporary variables as we update each time step
		tmpS = S[:]
		tmpA = A[:]

		# Divide R by S, update S = R - QS and A = QA + B
		[Q, r] = polyDiv(R, S)
		S = polySub(R, polyMult(Q, S))
		A = polyAdd(polyMult(Q, A), B)

		# Update R and B
		R = tmpS
		B = tmpA
		i += 1

	# Scale factor to make constant value of error locator be 1
	scale = A[0]

	# If scale becomes 0 at some point, divide by zero error and return None's
	if scale == 0:
		return [None, None]

	# Compute error locator by scaling
	errorLoc = [x / scale for x in A]

	# Compute error evaluator as S_initial * error locator (mod x^(2t))
	errorEval = polyMult(S0, errorLoc)
	[q, errorEval] = polyDiv(errorEval, GFPoly([0] * (2*t) + [1]))

	return [errorLoc, errorEval]

# Perform chien search - essentially check all elements in the field to find
# roots of the error locator polynomial
def chienSearch(errorLoc):
	results = []
	gamma = errorLoc
	for i in range(q-1):
		check = sum(gamma) % q
		if check == 0:
			results.append(GFInt(alpha).exponentiate(i))
		gamma = [gamma[j] * GFInt(alpha).exponentiate(j) \
					for j in range(len(gamma))]
	return results

# Take the formal derivative of a given polynomial using the power rule
def formalDeriv(poly):
	deriv = []
	for i in range(1, len(poly)):
		deriv.append(GFInt(i) * poly[i])
	return deriv

# Perform Forney's algorithm to actually calculate error locations and values
def forneyAlgo(errorEval, errorLoc, roots):
	# Locations are simply the inverses of the roots of the error locator
	locations = [logs[x.inverse()] for x in roots]

	# Evaluation points, take formal derivative of error locator for later
	evals = []
	deriv = formalDeriv(errorLoc)
	for r in roots:

		# Calculate error evaluator at roots
		num = [errorEval[j] * r.exponentiate(j) \
				for j in range(len(errorEval))]
		num = GFInt(sum(num) % q)

		# Calculate derivative of error locator at roots
		denom = [deriv[j] * r.exponentiate(j) for j in range(len(deriv))]
		denom = GFInt(sum(denom) % q)

		# If the derivative of the error locator evaluated at a root ends up
		# being 0, declare decoder failure and return None's to signify this
		if denom == 0:
			return [None, None]

		# Error value is -1 * error evalautor at root divided by derivative of
		# error locator at root
		evals.append(GFInt(-1) * num * denom.inverse())

	return [locations, evals]

# Decode y to a codeword
def decode(y):
	# Calculate syndrome
	syn = syndrome(y)

	# Calculate error locator and error evaluator polynomials
	[errorLoc, errorEval] = computeErrors(syn)

	# If either is None, we have a decoder failure, so return
	if errorLoc == None:
		return [None, None, y, True]

	# Find roots of error locator by performing a Chien Search
	roots = chienSearch(errorLoc)

	# If we don't have as many roots as the degree of the error locator, we
	# have a decoder failure
	if degree(errorLoc) != len(roots):
		fail = True

	# Calculate error locations and values
	[locations, evals] = forneyAlgo(errorEval, errorLoc, roots)

	# If locations is None, we have a decoder failure that is irrecoverable
	if locations == None:
		return [None, None, None, True]

	# Fix received vector y with locations and values
	fixed = y[:]
	fail = False
	for i in range(len(locations)):
		# If somehow we get a location out of bounds of y, decoder failure
		if locations[i] > len(y)-1:
			fail = True

		# Correct error at loc
		loc = locations[i]
		fixed[loc] -= evals[i]
	return [locations, evals, fixed, fail]

# Compute the generator polynomial
gener = computeGenerator()

# Compute SER and WER from a given message and probability for the channel
def computeErrorRates(m, prob):
	# Encode, go through the channel, and decode
	c = encode(m)
	y = channel(c, prob)
	[locations, evals, fixed, fail] = decode(y)

	# Default we have a word error
	wer = 1

	# If total decoder failure, all of the message bits are wrong as we cannot
	# decode at all and there is a word error
	if evals == None:
		return [k, 1]

	# If decoding succeded, no errors
	if not fail and fixed == c:
		return [0, 0]

	# SER is just the number of places where the initial message and decoded
	# message differ
	fixedM = fixed[-len(m):]
	ser = sum([int(fixedM[i] != m[i]) for i in range(len(m))])

	return [ser, wer]

# Generate a random message of length k
def genRandomMessage(k):
	m = []
	for i in range(k):
		m.append(random.randint(0, q-1))
	return GFPoly(m)

# Simulate running the code, probs and num_iter are lists of same length.
# probs[i] denotes a probability to run at and num_iter[i] denotes how many
# times to sample at that probability
def simulateCode(probs, num_iter, outfile):
	global p, m, q, t, n, k, alpha

	# File we write to
	out = open(outfile, "w")

	# Keep track of SER and WER
	SER = [0.0] * len(probs)
	WER = [0.0] * len(probs)

	# Simulate
	for i in range(len(num_iter)):
		prob = probs[i]
		n = num_iter[i]
		for j in range(n):
			m = genRandomMessage(k)
			[ser, wer] = computeErrorRates(m, prob)
			SER[i] += ser
			WER[i] += wer
		# Scale SER and WER as from computeErrorRates, they are just counts
		out.write("%f %f %f\n" % (prob, SER[i] / (n * k), WER[i] / n))
	return [SER, WER]

# Probabilities we are considering
probs = [0.01]
while probs[-1] < 1:
	probs.append(probs[-1] + 0.02)
probs = probs[:-1]

# How many times to run them
n_iter = [20000] * len(probs)

# Simulate
[SER, WER] = simulateCode(probs, n_iter, '47.txt')