#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <gmp.h>
#include <limits.h>

#define IS_PROBABLY_PRIME 1
#define IS_COMPOSITE 0
#define MILLER_RABBIN_K_PARAM 40

void sd_extract(mpz_t* s, mpz_t* d, mpz_t n);
void generate_key_components(mpz_t* rsa_modulus, mpz_t* public_exp, mpz_t* private_exp);
int miller_rabbin(mpz_t n, int k);

struct Public_Key {
	mpz_t pub_exp;
	mpz_t rsa_mod;
};

struct Private_Key {
	mpz_t priv_exp;
};

int main()
{
	struct Public_Key Pub_Key;
	struct Private_Key Priv_Key;
	mpz_t rsa_mod, pub_exp, priv_exp;

	mpz_init(rsa_mod);
	mpz_init(pub_exp);
	mpz_init(priv_exp);
	generate_key_components(&rsa_mod, &pub_exp, &priv_exp);


	mpz_init(Pub_Key.pub_exp);
	mpz_init(Pub_Key.rsa_mod);
	mpz_init(Priv_Key.priv_exp);

	mpz_set(Pub_Key.pub_exp, pub_exp);
	mpz_set(Pub_Key.rsa_mod, rsa_mod);

	mpz_set(Priv_Key.priv_exp, priv_exp);

	gmp_printf("Public Key public exponent: %Zd\n", Pub_Key.pub_exp);
	gmp_printf("Public Key RSA modulus: %Zd\n", Pub_Key.rsa_mod);
	gmp_printf("Public Key private exponent: %Zd\n", Priv_Key.priv_exp);
}


/* Generates the components for the public and private keys
 * in the RSA cryptosystem.
 *
 * Pass all arguments after mpz_init call.
  Arguments:
 * mpz_t* rsa_mod:  The RSA modulus which is used along with the
 *                  public and private exponents to form the public
 *                  and private keys.
 * mpz_t* pub_exp:  The public exponent that forms the public key
 *                  along with the rsa_modulus.
 * mpz_t* priv_exp: The private exponent that form the private key
 *                  along with the rsa_modulus.
 *
 * Returns:
 * Void (Sets all pointers at the end of the function).
 */
void generate_key_components(mpz_t* rsa_mod, mpz_t* pub_exp, mpz_t* priv_exp)
{
	/* Randomization */
	gmp_randstate_t rand_state; /* Randomized state used to generate p and q. */
	int seed  = time(NULL);

	/*
	 * Used in prime generation.
	 *
	 * p, q: The prime factors of rsa_mod.
	 * p_max: The max value of p.
	 * q_max: The max value, (to-be) offset in size by p_q_offset
	 * TO BE IMPLEMENTED: p_q_offset: Randomized integer, used to offset the max value of
	 *             q from p to make factorization harder.
	 */
	mpz_t p, q, p_max, q_max;
	int p_is_prime = 0, q_is_prime = 0;
	/* int p_q_offset; */

	/* Variables for computing and storing the Carmichael function of RSA modulus. */
	mpz_t carmichael_function, p_sub_1, q_sub_1;

	/* Variables used to calculate the public (or "encryption") exponent. */
	mpz_t temp_pub_exp, pub_exp_gcd, carmichael_function_sub_2;
	int pub_exp_found = 0;
	
	/* Initialize all prime generation variables. */
	mpz_init(p);
	mpz_init(q);
	mpz_init(p_max);
	mpz_init(q_max);
	mpz_set_ui(p_max, 10000);
	mpz_set_ui(q_max, 10000);

	/* Initiates variables for computing Carmichael function. */
	mpz_init(carmichael_function);
	mpz_init(p_sub_1);
	mpz_init(q_sub_1);

	/* Initiates variables for computing the public key exponent. */
	mpz_init(temp_pub_exp);
	mpz_init(pub_exp_gcd);
	mpz_init(carmichael_function_sub_2);
	
	/* Setting up random state and seed. */
	gmp_randinit_mt(rand_state);
	gmp_randseed_ui(rand_state, seed);

	/* Generate the prime factors for the RSA modulus.
	 * NOTE: At the moment we generate only up to max unsigned int for testing
	 * because this implementation of Miller-Rabbin seems a bit slow at the moment.
	 */
	while (!p_is_prime)
	{
		mpz_urandomm(p, rand_state, p_max);
		if (mpz_odd_p(p) != 0)
		{
			gmp_printf("Randomized p: %Zd\n", p);
			p_is_prime = miller_rabbin(p, MILLER_RABBIN_K_PARAM);
		}
	}

	while (!q_is_prime)
	{
		mpz_urandomm(q, rand_state, p_max);
		if (mpz_odd_p(q) != 0)
		{
			gmp_printf("Randomized q: %Zd\n", q);
			q_is_prime = miller_rabbin(p, MILLER_RABBIN_K_PARAM);
		}
	}

	gmp_printf("p: %Zd\nq: %Zd\n", p, q);

	/* Compute the RSA modulus. */
	mpz_mul(*rsa_mod, p, q);
	gmp_printf("The rsa modulus is: %Zd\n", rsa_mod);

	/* Calculate Carmichael's totient function. */
	mpz_sub_ui(p_sub_1, p, 1);
	mpz_sub_ui(q_sub_1, q, 1);
	mpz_lcm(carmichael_function, p_sub_1, q_sub_1);
	gmp_printf("The carmichael function  is: %Zd\n", carmichael_function);

	/* Calculate the public key exponent.
	 * Randomly select a number in (1, carmichael_function), check if
	 * the number and carmichael_function are coprime, if so , set the
	 * public exponent to be that number.
	 */
	while (!pub_exp_found)
	{
		mpz_sub_ui(carmichael_function_sub_2, carmichael_function, 2);
		mpz_urandomm(temp_pub_exp, rand_state, carmichael_function_sub_2);
		mpz_add_ui(temp_pub_exp, temp_pub_exp, 2);
		mpz_gcd(pub_exp_gcd, temp_pub_exp, carmichael_function);
		pub_exp_found = (mpz_cmp_ui(pub_exp_gcd, 1) ==  0);
	}

	mpz_set(*pub_exp, temp_pub_exp);
	gmp_printf("The public key exponent: %Zd\n", *pub_exp);

	/* Calculate the private key exponent, which is the modular multiplicative inverse
	 * of pub_exp modulo carmichael_function. */
	mpz_invert(*priv_exp, *pub_exp, carmichael_function);
	gmp_printf("The private key exponent: %Zd\n", *priv_exp);
}

/* The extended Euclidean algorithm is used to calculate the private key exponent,
 * which is the modular multiplicative inverse of "pub_exp modulo carmichael_function".
 *
 * Arguments:
 * mpz_t* priv_exp: A pointer to the variable which is our private key exponent.
 * mpz_t pub_exp: The public key exponent.
 * mpz_t carmichael_function : The carmichael function.
 *
 * Returns:
 * void */
/*
void extended_euclidean(mpz_t* priv_exp, mpz_t pub_exp, mpz_t carmichael_function)
{
	mpz_t temp_val, r, old_r, quotient, s, old_s, t, old_t;

	mpz_init(temp_val);
	mpz_init(r);
	mpz_init(old_r);
	mpz_init(quotient);
	mpz_init(s);
	mpz_init(old_s);
	mpz_init(t);
	mpz_init(old_t);

	mpz_set(old_r, a);
	mpz_set(r, b);
	mpz_set_ui(old_s, 1);
	mpz_set_ui(s, 0);
	mpz_set_ui(old_t, 0);
	mpz_set_ui(t, 1);
	
	while (true)
	{
	
		while(mpz_compare_ui(r,0) == 0)
		{
			quotient = mpz_div();
		}
	}

}
*/
/* Tells us whether a number is prime using the probabilistic
 * version of the Miller-Rabbin primality test.
 *
 * Arguments:
 * mpz_t n: An odd value > 3 we are performing the primality test on.
 * mpz_t k: A parameter that determines how accurate the test
 *        will be.
 *
 * Returns:
 * 1: If n is "probably prime"
 * 0: If n is a composite number.
 */
int miller_rabbin(mpz_t n, int k)
{
	/* Iter. vars. */
	int i;
	mpz_t j;

	/* Booleans */
	int continue_witness = 0;
	int mpz_x_eq_1;
	int mpz_x_eq_n_sub_1;

	/* "Helper vars". */
	mpz_t n_sub_1, n_sub_3, s_sub_1;

	/* Major parameters of the Miller-Rabbin algorithm. */
	mpz_t x, a, s, d;

	gmp_randstate_t rand_state; /* Randomized state used to generate a. */

	/* Initiate all gmp variables to zero. */
	mpz_init(j);
	mpz_init(n_sub_1);
	mpz_init(n_sub_3);
	mpz_init(s_sub_1);
	mpz_init(x);
	mpz_init(a);
	mpz_init(s);
	mpz_init(d);

	mpz_sub_ui(n_sub_1, n, 1);
	mpz_sub_ui(n_sub_3, n, 3);

	gmp_randinit_mt(rand_state);

	/* Extract the s and d parameters. */
	sd_extract(&s, &d, n);

	mpz_sub_ui(s_sub_1, s, 1);

	for (i = 0; i < k; i++)
	{
		continue_witness = 0;

		/* Initialize the random integer a. */
		mpz_urandomm(a, rand_state, n_sub_3);
		mpz_add_ui(a, a, 2);

		/* Initialize x. */
		mpz_pow_ui(x, a, mpz_get_ui(d));
		mpz_mod(x, x, n);

		mpz_x_eq_1 = (mpz_cmp_ui(x,1) == 0);
		mpz_x_eq_n_sub_1 = (mpz_cmp(x, n_sub_1) == 0);
		if (mpz_x_eq_1 || mpz_x_eq_n_sub_1)
		{
			continue;
		}
		for (mpz_init(j); mpz_cmp(j, s_sub_1) < 0; mpz_add_ui(j, j, 1))
		{
			mpz_pow_ui(x, x, 2);
			mpz_mod(x, x, n);
			mpz_x_eq_1 = (mpz_cmp_ui(x,1) == 0);
			mpz_x_eq_n_sub_1 = (mpz_cmp(x, n_sub_1) == 0);

			if (mpz_x_eq_1)
				return IS_COMPOSITE;
			else if (mpz_x_eq_n_sub_1)
			{
				continue_witness = 1;
				break;
			}
		}
		if (continue_witness)
			continue;
		return IS_COMPOSITE;
	}
	return IS_PROBABLY_PRIME;
}

/* Sets positive integer pointers s and d, d odd, where:
 * (2^s * d) = (n-1).
 * This is used for the Miller-Rabbin primality test. */
void sd_extract(mpz_t* s, mpz_t* d, mpz_t n)
{
	int found_sd = 0, n_diff_too_big;

	mpz_t n_diff; /* Simply n-1 in a separate value. */
	mpz_t temp_n_diff;
	mpz_t temp_n_diff_base;
	mpz_t temp_d, temp_s;

	/* Initiate all gmp variables. */
	mpz_init(n_diff);
	mpz_init(temp_n_diff);
	mpz_init(temp_n_diff_base);
	mpz_init(temp_s);
	mpz_init(temp_d);

	/* Set initial values of a few variables. */
	mpz_sub_ui(n_diff, n, 1);
	mpz_set_ui(temp_n_diff_base, 2);
	mpz_set_ui(temp_d, 1);

	while (!found_sd)
	{
		mpz_set_ui(temp_s, 1);
		n_diff_too_big = 0;
		while ((!found_sd) && (!n_diff_too_big))
		{
			mpz_pow_ui(temp_n_diff, temp_n_diff_base, mpz_get_ui(temp_s));
			mpz_mul(temp_n_diff, temp_n_diff, temp_d);
			found_sd = (mpz_cmp(temp_n_diff, n_diff) == 0);
			n_diff_too_big = (mpz_cmp(temp_n_diff, n_diff) > 0);
			if (!found_sd)
				mpz_add_ui(temp_s, temp_s, 1);
		}
		if (!found_sd)
			mpz_add_ui(temp_d, temp_d, 1);
	}
	mpz_set(*d, temp_d);
	mpz_set(*s, temp_s);
}
