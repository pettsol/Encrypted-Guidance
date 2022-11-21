#ifndef HE_GUIDANCE_H
#define HE_GUIDANCE_H

#include <gmp.h>
#include <gmpxx.h>
#include <stdint.h>
#include <math.h>
#include <queue>

class Encrypted_ilos_guidance
{
	private:
		std::queue<mpz_class> encrypted_waypoints;
		std::queue<mpz_class> courses;
		std::queue<mpz_class> encrypted_pi;
		mpz_t g;
		mpz_t N;
		mpz_t N2;
		mpz_t lambda;
		mpz_t p2;
		mpz_t q2;
		mpz_t mu;
		mpz_t gamma_p;
		mpz_t gamma_kp;
		mpz_t gamma_ki;
		mpz_t gamma_inv_trig;
		mpz_t gamma_inverse;
		mpz_t gamma_time;
		mpz_t kp_bar;
		mpz_t ki_bar;
		mpz_t delta_t;
		mpz_t c_integral;
		mpz_t current_origin_x;
		mpz_t current_origin_y;
		mpz_t next_wp_x;
		mpz_t next_wp_y;
		mpz_t current_sin_pi;
		mpz_t current_cos_pi;
		mpz_t current_c_pi;
		gmp_randstate_t randstate;
		uint32_t msgsize;
		float threshold;
	public:
		Encrypted_ilos_guidance(
				const mpz_t N_in,
				const mpz_t N2_in,
				const mpz_t lambda_in,
				const mpz_t g_in,
				const mpz_t p2_in,
				const mpz_t q2_in,
				const mpz_t mu,
				const mpz_t gamma_p_in,
				const mpz_t gamma_kp_in,
				const mpz_t gamma_ki_in,
				const mpz_t gamma_inverse_in,
				const mpz_t gamma_inv_trig_in,
				const mpz_t gamma_time_in,
				const float threshold,
				const float kp,
				const float ki,
				const float delta_t_in,
				const uint32_t size);
		void preprocessing(
				float *waypoints,
				size_t n); 
		void quantize_and_encrypt(
				mpz_t c_x,
				mpz_t c_y,
				const float x_n,
				const float y_n);
		void iterate(
				mpz_t c_psi_d,
                                mpz_t c_xe,
                                mpz_t c_ye,
				uint32_t *b,
                                const mpz_t c_xn,
                                const mpz_t c_yn); // Iterate ILOS guidance law
		void decrypt_and_recover(
				mpf_t psi_d_f,
				uint32_t *b,
				const mpz_t c_psi_d,
				const mpz_t c_x_e,
				const mpz_t c_y_e);
		void rho(
				mpz_t out,
				const mpf_t in,
				const mpz_t gamma,
				const mpz_t ptspace);
		void rho_inv(
				mpf_t out,
				const mpz_t in,
				const mpz_t gamma,
				const mpz_t ptspace);
};
#endif
