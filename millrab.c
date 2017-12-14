#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <gmp.h>

#define IS_PROBABLY_PRIME 1
#define IS_COMPOSITE 0
#define MILLER_RABBIN_K_PARAM 40

void sd_extract(mpz_t* s, mpz_t* d, mpz_t n);
int miller_rabbin(mpz_t n, int k);

int main()
{
	int miller_rabbin_result;
	
	mpz_t n; /* The prime we want to test. */
	
	/* Two significant parameters for Miller-Rabbin. */
	mpz_t s;
	mpz_t d;

	int k = 40; /* A parameter deciding accuracy for Miller-Rabbin. */
	
	/* Initiate all gmp variables. */
	mpz_init(n);
	mpz_init(s);
	mpz_init(d);

	mpz_set_str(n, "15486489", 10);
	
	gmp_printf("n: %Zd\n", n);

	sd_extract(&s, &d, n);

	gmp_printf("s: %Zd and d: %Zd\n", s, d);

	miller_rabbin_result = miller_rabbin(n, k);

	printf("Result Miller-Rabbin: %d\n", miller_rabbin_result);
}


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
 * */
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
