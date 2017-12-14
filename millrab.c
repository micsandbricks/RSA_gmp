#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <gmp.h>

#define IS_PROBABLY_PRIME 1
#define IS_COMPOSITE 0
#define MILLER_RABBIN_K_PARAM 40

void sd_extract(mpz_t* s, mpz_t* d, mpz_t n);

int main()
{
	mpz_t n;
	mpz_t s;
	mpz_t d;

	/* Initiate all gmp variables. */
	mpz_init(n);
	mpz_init(s);
	mpz_init(d);

	mpz_set_str(n, "179424793", 10);
	
	gmp_printf("n: %Zd\n", n);

	sd_extract(&s, &d, n);

	gmp_printf("s: %Zd and d: %Zd\n", s, d);
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
