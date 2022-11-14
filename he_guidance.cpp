#include <algorithm>
#include <iostream>

#include "he_guidance.h"
#include "joye_libert_journal/joye_libert.h"

Encrypted_ilos_guidance::Encrypted_ilos_guidance(
						const mpz_t N_in,
						const mpz_t y_in,
						const mpz_t p_in,
						const mpz_t gamma_p_in,
						const mpz_t gamma_kp_in,
						const mpz_t gamma_ki_in,
						const mpz_t gamma_inverse_in,
						const mpz_t gamma_inv_trig_in,
						const mpz_t gamma_time_in,
						const float threshold_in,
						const float kp,
						const float ki,
						const float delta_t_in,
						const uint32_t size)
{
	// Constructor;
	// Init all stuff
	mpz_init_set(N, N_in);
	mpz_init_set(y, y_in);
	mpz_init_set(p, p_in);
	mpz_init_set(gamma_p, gamma_p_in);
	mpz_init_set(gamma_kp, gamma_kp_in);
	mpz_init_set(gamma_ki, gamma_ki_in);
	mpz_init_set(gamma_inverse, gamma_inverse_in);
	mpz_init_set(gamma_inv_trig, gamma_inv_trig_in);
	mpz_init_set(gamma_time, gamma_time_in);
	mpz_init(kp_bar);
	mpz_init(ki_bar);
	mpz_init(delta_t);
	mpz_init(c_integral);
	mpz_init(next_wp_x);
	mpz_init(next_wp_y);
	mpz_init(current_origin_x);
	mpz_init(current_origin_y);
	mpz_init(current_sin_pi);
	mpz_init(current_cos_pi);
	mpz_init(current_c_pi);
	gmp_randinit_mt(randstate);
	msgsize = size;
	threshold = threshold_in;

	mpf_t kp_f, ki_f, delta_t_f;
	mpf_init_set_d(kp_f, kp);
	mpf_init_set_d(ki_f, ki);
	mpf_init_set_d(delta_t_f, delta_t_in);
	rho(kp_bar, kp_f, gamma_kp);
	rho(ki_bar, ki_f, gamma_ki);
	rho(delta_t, delta_t_f, gamma_time);
}

void Encrypted_ilos_guidance::preprocessing(float *waypoints, size_t n)
{
	float x_prev, y_prev, x_next, y_next, pi;
	mpf_t x_f, y_f, sin_pi, cos_pi, pi_f;
	mpf_init(x_f);
	mpf_init(y_f);
	mpf_init(sin_pi);
	mpf_init(cos_pi);
	mpf_init(pi_f);
	mpz_t c_x, c_y, x_bar, y_bar, sin_pi_bar, cos_pi_bar, c_pi;
	mpz_init(c_x);
	mpz_init(c_y);
	mpz_init(x_bar);
	mpz_init(y_bar);
	mpz_init(sin_pi_bar);
	mpz_init(cos_pi_bar);
	mpz_init(c_pi);
	
	for (int i = 1; i < n; i++)
	{
		x_prev = *(waypoints++);
		y_prev = *(waypoints++);
		x_next = *(waypoints);
		y_next = *(waypoints+1);

		pi = atan2(y_next - y_prev, x_next - x_prev);

		// NOTE: The second waypoint is the "first waypoint" we 
		// 	 are going to visit
		mpf_set_d(x_f, x_prev);
		mpf_set_d(y_f, y_prev);
		mpf_set_d(sin_pi, sin(pi));
		mpf_set_d(cos_pi, cos(pi));

		rho(x_bar, x_f, gamma_p);
		rho(y_bar, y_f, gamma_p);
		rho(sin_pi_bar, sin_pi, gamma_p);
		rho(cos_pi_bar, cos_pi, gamma_p);

		joye_libert_encrypt(c_x, randstate, x_bar, y, N, msgsize);
		joye_libert_encrypt(c_y, randstate, y_bar, y, N, msgsize);

		mpz_class c_x_class(c_x), c_y_class(c_y), sin_pi_class(sin_pi_bar), cos_pi_class(cos_pi_bar);
		encrypted_waypoints.push(c_x_class);
		encrypted_waypoints.push(c_y_class);
		courses.push(sin_pi_class);
		courses.push(cos_pi_class);

		mpf_set_d(pi_f, pi);
		rho(c_pi, pi_f, gamma_inverse);
		joye_libert_encrypt(c_pi, randstate, c_pi, y, N, msgsize);
		mpz_class c_pi_class(c_pi);
		encrypted_pi.push(c_pi_class);


	}
	// Add the last waypoint
	x_next = *(waypoints++);
	y_next = *(waypoints++);
	mpf_set_d(x_f, x_next);
	mpf_set_d(y_f, y_next);
	rho(x_bar, x_f, gamma_p);
	rho(y_bar, y_f, gamma_p);

	joye_libert_encrypt(c_x, randstate, x_bar, y, N, msgsize);

	joye_libert_encrypt(c_y, randstate, y_bar, y, N, msgsize);

	mpz_class c_x_class(c_x), c_y_class(c_y);
	encrypted_waypoints.push(c_x_class);
	encrypted_waypoints.push(c_y_class);
	
	// Extract first waypoint as origin, and second waypoint and next WP
	mpz_class current_origin_x_class = encrypted_waypoints.front(); encrypted_waypoints.pop();
	mpz_class current_origin_y_class = encrypted_waypoints.front(); encrypted_waypoints.pop();
	mpz_class next_wp_x_class = encrypted_waypoints.front(); encrypted_waypoints.pop();
	mpz_class next_wp_y_class = encrypted_waypoints.front(); encrypted_waypoints.pop();
	mpz_class current_sin_pi_class = courses.front(); courses.pop();
	mpz_class current_cos_pi_class = courses.front(); courses.pop();
	mpz_class current_c_pi_class = encrypted_pi.front(); encrypted_pi.pop();

	mpz_set(current_origin_x, current_origin_x_class.get_mpz_t());
	mpz_set(current_origin_y, current_origin_y_class.get_mpz_t());
	mpz_set(next_wp_x, next_wp_x_class.get_mpz_t());
	mpz_set(next_wp_y, next_wp_y_class.get_mpz_t());
	mpz_set(current_sin_pi, current_sin_pi_class.get_mpz_t());
	mpz_set(current_cos_pi, current_cos_pi_class.get_mpz_t());
	mpz_set(current_c_pi, current_c_pi_class.get_mpz_t());

	mpz_t c_int_0;
	mpz_init_set_ui(c_int_0, 0);
	joye_libert_encrypt(c_int_0, randstate, c_int_0, y, N, msgsize);
	
	mpz_set(c_integral, c_int_0);
}

void Encrypted_ilos_guidance::quantize_and_encrypt(
						mpz_t c_x,
						mpz_t c_y,
						const float x_n,
						const float y_n)
{
	mpz_t x_n_bar, y_n_bar;

	mpz_init(x_n_bar);
	mpz_init(y_n_bar);

	mpf_t x_n_f, y_n_f;
	mpf_init_set_d(x_n_f, x_n);
	mpf_init_set_d(y_n_f, y_n);

	rho(x_n_bar, x_n_f, gamma_p);
	rho(y_n_bar, y_n_f, gamma_p);

	mpz_init(c_x);
	mpz_init(c_y);

	joye_libert_encrypt(c_x, randstate, x_n_bar, y, N, msgsize);
	joye_libert_encrypt(c_y, randstate, y_n_bar, y, N, msgsize);
}

void Encrypted_ilos_guidance::iterate(
				mpz_t c_psi_d,
			       	mpz_t c_xe,
			       	mpz_t c_ye,
			       	uint32_t *b,
				const mpz_t c_xn,
			       	const mpz_t c_yn)
{
	
	// If waypoint reached, get next waypoint or terminate
	if (*b)
	{
		
		std::cout << "Waypoint reached!\n";
		(*b) = 0;
		if (encrypted_waypoints.empty())
		{
			std::cout << "Final destination reached!\n";
			return;
		}
		
		// Go over this again - Clearly we should start at the second waypoint
		mpz_class next_wp_x_class = encrypted_waypoints.front(); encrypted_waypoints.pop();
		mpz_class next_wp_y_class = encrypted_waypoints.front(); encrypted_waypoints.pop();
		mpz_class next_sin_pi_class = courses.front(); courses.pop();
		mpz_class next_cos_pi_class = courses.front(); courses.pop();
		mpz_class next_c_pi_class = encrypted_pi.front(); encrypted_pi.pop();

		mpz_set(current_origin_x, next_wp_x);
		mpz_set(current_origin_y, next_wp_y);
		mpz_set(next_wp_x, next_wp_x_class.get_mpz_t());
		mpz_set(next_wp_y, next_wp_y_class.get_mpz_t());
		mpz_set(current_sin_pi, next_sin_pi_class.get_mpz_t());
		mpz_set(current_cos_pi, next_cos_pi_class.get_mpz_t());
		mpz_set(current_c_pi, next_c_pi_class.get_mpz_t());

	}

	// Compute encrypted ilos guidance
	mpz_t tmp_x, tmp_y, inv, c_yep, c_xep;
	mpz_init(tmp_x);
	mpz_init(tmp_y);
	mpz_init(inv);
	mpz_init(c_yep);
	mpz_init(c_xep);

	// Compute encrypted cross-track error
	// homomorphically
	mpz_invert(inv, current_origin_x, N);
	mpz_mul(c_yep, c_xn, inv);
	mpz_mod(c_yep, c_yep, N);

	mpz_powm(c_yep, c_yep, current_sin_pi, N);
	mpz_invert(c_yep, c_yep, N);

	mpz_invert(inv, current_origin_y, N);
	mpz_mul(tmp_x, c_yn, inv);
	mpz_mod(tmp_x, tmp_x, N);
	mpz_powm(tmp_x, tmp_x, current_cos_pi, N);

	mpz_mul(c_yep, c_yep, tmp_x);
	mpz_mod(c_yep, c_yep, N);
	
	// Update integral
	mpz_powm(tmp_x, c_yep, delta_t, N);
	mpz_mul(c_integral, c_integral, tmp_x);
	mpz_mod(c_integral, c_integral, N);

	// Compute c_psi_d
	mpz_powm(tmp_x, c_yep, kp_bar, N);
	mpz_invert(tmp_x, tmp_x, N);
	mpz_powm(inv, c_integral, ki_bar, N);
	mpz_invert(inv, inv, N);
	mpz_mul(c_psi_d, current_c_pi, tmp_x);
	mpz_mod(c_psi_d, c_psi_d, N);
	mpz_mul(c_psi_d, c_psi_d, inv);
	mpz_mod(c_psi_d, c_psi_d, N);

	// Compute encrypted along-track and cross-track distance to the next waypoint
	// First part
	mpz_t x_tilde;
	mpz_init(x_tilde);
	//mpz_invert(inv, c_xn, N);
	mpz_invert(inv, next_wp_x, N);
	//mpz_mul(x_tilde, next_wp_x, inv);
	mpz_mul(x_tilde, c_xn, inv);
	mpz_mod(x_tilde, x_tilde, N);
	mpz_powm(c_yep, x_tilde, current_sin_pi, N);
	mpz_powm(c_xep, x_tilde, current_cos_pi, N);
	// y-coordinate is to be subtracted, so we need the multiplicative
	// inverse
	mpz_invert(c_yep, c_yep, N);

	// Second part
	mpz_invert(inv, next_wp_y, N);
	mpz_mul(tmp_y, c_yn, inv);
	//mpz_invert(inv, c_yn, N);
	//mpz_mul(tmp_y, next_wp_y, inv);
	mpz_mod(tmp_y, tmp_y, N);
	mpz_powm(tmp_x, tmp_y, current_sin_pi, N);
	mpz_powm(tmp_y, tmp_y, current_cos_pi, N);

	// Homomorphically add first and second part
	mpz_mul(c_xe, c_xep, tmp_x);
	mpz_mul(c_ye, c_yep, tmp_y);
	mpz_mod(c_xe, c_xe, N);
	mpz_mod(c_ye, c_ye, N);
}

void Encrypted_ilos_guidance::decrypt_and_recover(
						mpf_t psi_d_f,
						uint32_t *b,
						const mpz_t c_psi_d,
					        const mpz_t c_x_e,
						const mpz_t c_y_e)
{
	(*b) = 0;

	mpz_t psi_d_bar, x_e_bar, y_e_bar;
	mpz_init(psi_d_bar);
	mpz_init(x_e_bar);
	mpz_init(y_e_bar);

	joye_libert_decrypt(psi_d_bar, c_psi_d, p, y, msgsize);
	joye_libert_decrypt(x_e_bar, c_x_e, p, y, msgsize);
	joye_libert_decrypt(y_e_bar, c_y_e, p, y, msgsize);

	mpf_t x_e_p, y_e_p;
	mpf_init(x_e_p);
	mpf_init(y_e_p);

	// Set the right gamma!	
	rho_inv(psi_d_f, psi_d_bar, gamma_inverse);
	rho_inv(x_e_p, x_e_bar, gamma_inv_trig);
	rho_inv(y_e_p, y_e_bar, gamma_inv_trig);

	float x_p_abs = abs(mpf_get_d(x_e_p));
	float y_p_abs = abs(mpf_get_d(y_e_p));

	float delta = std::max(x_p_abs, y_p_abs);

	gmp_printf("x_e_p: %Ff\n", x_e_p);
	gmp_printf("y_e_p: %Ff\n", y_e_p);

	if (delta < threshold)
	{
		std::cout << "Delta less than threshold, setting b = 1\n";
		(*b) = 1;
	}
	// Update the shared variable
	// Send phi_d to the control system?
}

void Encrypted_ilos_guidance::rho(mpz_t out, mpf_t in, mpz_t gamma)
{
	mpf_t tmp;
	mpf_init(tmp);
	mpf_set_z(tmp, gamma);
	mpf_mul(tmp, in, tmp);

	mpz_t size;
	mpz_init(size);
	mpz_ui_pow_ui(size, 2, msgsize);

	mpz_set_f(out, tmp);
	mpz_mod(out, out, size);
}

void Encrypted_ilos_guidance::rho_inv(mpf_t out, mpz_t in, mpz_t gamma)
{
	
	mpz_t size, halfsize;
	mpz_init(size);
	mpz_init(halfsize);
	mpz_ui_pow_ui(size, 2, msgsize);
	mpz_ui_pow_ui(halfsize, 2, msgsize-1);

	
	mpz_t test; 
	mpz_init(test);
	
	mpz_sub(test, in, halfsize);

	if (mpz_sgn(test) != -1)
	{
		// negative number
		mpz_sub(in, in, size);
	}

	mpf_t gamma_f;
	mpf_init(gamma_f);

	mpf_set_z(gamma_f, gamma);
	mpf_set_z(out, in);
	mpf_div(out, out, gamma_f);
}

