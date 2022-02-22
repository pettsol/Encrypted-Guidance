#include <algorithm>
#include <iostream>

#include "he_guidance.h"
#include "joye_libert_journal/joye_libert.h"

Encrypted_ilos_guidance::Encrypted_ilos_guidance(
						const mpz_t N_in,
						const mpz_t y_in,
						const mpz_t p_in,
						const mpz_t gamma_in,
						const mpz_t gamma_inverse_in,
						const mpz_t gamma_inv_trig_in,
						const mpz_t gamma_time_in,
						const float threshold_in,
						const float kp,
						const float ki,
						const float delta_t_in,
						const uint32_t size)
{
	//?? Constructor;
	// Init all stuff
	//
	// TODO: TA HØYDE FOR AT REGULATORPARAMETRE MÅ SKALERES
	mpz_init_set(N, N_in);
	mpz_init_set(y, y_in);
	mpz_init_set(p, p_in);
	mpz_init_set(gamma, gamma_in);
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
	rho(kp_bar, kp_f, gamma);
	rho(ki_bar, ki_f, gamma);
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

		std::cout << "pi: " << pi << std::endl;
		// NOTE: The second waypoint is the "first waypoint" we 
		// 	 are going to visit
		mpf_set_d(x_f, x_prev);
		mpf_set_d(y_f, y_prev);
		mpf_set_d(sin_pi, sin(pi));
		mpf_set_d(cos_pi, cos(pi));

		rho(x_bar, x_f, gamma);
		rho(y_bar, y_f, gamma);
		rho(sin_pi_bar, sin_pi, gamma);
		rho(cos_pi_bar, cos_pi, gamma);

		gmp_printf("sin_pi_bar: %Zd\n", sin_pi_bar);
		gmp_printf("cos_pi_bar: %Zd\n", cos_pi_bar);

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
	std::cout << "Post loop\n";
	// Add the last waypoint
	x_next = *(waypoints++);
	y_next = *(waypoints++);
	std::cout << "xy = " << x_next << y_next << std::endl;
	mpf_set_d(x_f, x_next);
	mpf_set_d(y_f, y_next);
	rho(x_bar, x_f, gamma);
	rho(y_bar, y_f, gamma);

	gmp_printf("x_bar: %Zd\n", x_bar);

	joye_libert_encrypt(c_x, randstate, x_bar, y, N, msgsize);

	gmp_printf("c_x: %Zd\n", c_x);

	joye_libert_encrypt(c_y, randstate, y_bar, y, N, msgsize);
	// DEBUG
	gmp_printf("x_f: %Ff\n", x_f);
	joye_libert_decrypt(x_bar, c_x, p, y, msgsize);
	rho_inv(x_f, x_bar, gamma);
	gmp_printf("x_f: %Ff\n", x_f);
	// END DEBUG

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

	std::cout << "Before current origin\n";
	mpz_set(current_origin_x, current_origin_x_class.get_mpz_t());
	mpz_set(current_origin_y, current_origin_y_class.get_mpz_t());
	mpz_set(next_wp_x, next_wp_x_class.get_mpz_t());
	mpz_set(next_wp_y, next_wp_y_class.get_mpz_t());
	mpz_set(current_sin_pi, current_sin_pi_class.get_mpz_t());
	mpz_set(current_cos_pi, current_cos_pi_class.get_mpz_t());
	mpz_set(current_c_pi, current_c_pi_class.get_mpz_t());

	std::cout << "Before integral\n";
	mpz_t c_int_0;
	mpz_init_set_ui(c_int_0, 0);
	joye_libert_encrypt(c_int_0, randstate, c_int_0, y, N, msgsize);
	
	std::cout << "Before set c_integral\n";
	mpz_set(c_integral, c_int_0);

	std::cout << "After set c_integral\n";
}

void Encrypted_ilos_guidance::quantize_and_encrypt(
						mpz_t c_x,
						mpz_t c_y,
						const float x_n,
						const float y_n,
						uint32_t *b)
{
	mpz_t x_n_bar, y_n_bar;

	mpz_init(x_n_bar);
	mpz_init(y_n_bar);

	mpf_t x_n_f, y_n_f;
	mpf_init_set_d(x_n_f, x_n);
	mpf_init_set_d(y_n_f, y_n);

	rho(x_n_bar, x_n_f, gamma);
	rho(y_n_bar, y_n_f, gamma);

	mpz_init(c_x);
	mpz_init(c_y);

	joye_libert_encrypt(c_x, randstate, x_n_bar, y, N, msgsize);
	joye_libert_encrypt(c_y, randstate, y_n_bar, y, N, msgsize);
}

// Use a class for this? Queue of encrypted waypoints, angles, ... ?
void Encrypted_ilos_guidance::iterate(
				mpz_t c_psi_d,
			       	mpz_t c_xe,
			       	mpz_t c_ye,
			       	uint32_t *b,
				const mpz_t c_xn,
			       	const mpz_t c_yn)
{
	
	// If waypoint reached, get next waypoint or terminate
	std::cout << "b: " << *b << std::endl;
	if (*b)
	{
		
		std::cout << "Inside logic check for b\n";
		(*b) = 0;
		if (encrypted_waypoints.empty())
		{
			std::cout << "Inside empty\n";
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

		// DEBUG
		mpz_t x_debug, y_debug, pi_debug;
		mpz_init(x_debug);
		mpz_init(y_debug);
		mpz_init(pi_debug);

		mpf_t x_debug_f, y_debug_f, pi_debug_f;
		mpf_init(x_debug_f);
		mpf_init(y_debug_f);
		mpf_init(pi_debug_f);

		// DECRYPT AND DEBUG //
		joye_libert_decrypt(x_debug, current_origin_x, p, y, msgsize);
		joye_libert_decrypt(y_debug, current_origin_y, p, y, msgsize);
		joye_libert_decrypt(pi_debug, current_c_pi, p, y, msgsize);
		
		rho_inv(x_debug_f, x_debug, gamma);
		rho_inv(y_debug_f, y_debug, gamma);
		rho_inv(pi_debug_f, pi_debug, gamma);
		
		gmp_printf("Current origin x: %Ff\n", x_debug_f);
		gmp_printf("Current origin y: %Ff\n", y_debug_f);
		gmp_printf("Current pi: %Ff\n", pi_debug_f);
		
		joye_libert_decrypt(x_debug, next_wp_x, p, y, msgsize);
		joye_libert_decrypt(y_debug, next_wp_y, p, y, msgsize);
		
		rho_inv(x_debug_f, x_debug, gamma);
		rho_inv(y_debug_f, y_debug, gamma);
		
		gmp_printf("Next wp x: %Ff\n", x_debug_f);
		gmp_printf("Next wp y: %Ff\n", y_debug_f);

		//rho_inv(debug_1f, debug_1, gamma_inverse);
		//gmp_printf("Output cross-track error: %Ff\n", debug_1f);
	// END DEBUG

		// END DEBUG
	}
	// Coordinates begin in the PREVIOUS waypoint, right?
	// 	Establish variables for the previous waypoint,
	// 	i.e., the origin of the current path-fixed frame?

	// Compute encrypted ilos guidance
	mpz_t tmp_x, tmp_y, inv, c_yep, c_xep;
	mpz_init(tmp_x);
	mpz_init(tmp_y);
	mpz_init(inv);
	mpz_init(c_yep);
	mpz_init(c_xep);

	// Compute cross-track error
	mpz_invert(inv, current_origin_x, N);
	mpz_mul(c_yep, c_xn, inv);
	mpz_mod(c_yep, c_yep, N);

	// Sjekk at subtraksjon virker
	mpz_t debug_1;
	mpz_init(debug_1);
	mpf_t debug_1f;
	mpf_init(debug_1f);
	joye_libert_decrypt(debug_1, c_yep, p, y, msgsize);

	rho_inv(debug_1f, debug_1, gamma);
	gmp_printf("Output error in NED x-position: %Ff\n", debug_1f);

	// End check

	mpz_powm(c_yep, c_yep, current_sin_pi, N);
	mpz_invert(c_yep, c_yep, N);

	mpz_invert(inv, current_origin_y, N);
	mpz_mul(tmp_x, c_yn, inv);
	mpz_mod(tmp_x, tmp_x, N);
	mpz_powm(tmp_x, tmp_x, current_cos_pi, N);

	mpz_mul(c_yep, c_yep, tmp_x);
	mpz_mod(c_yep, c_yep, N);
	
	// DECRYPT AND DEBUG //
	joye_libert_decrypt(debug_1, c_yep, p, y, msgsize);
	
	rho_inv(debug_1f, debug_1, gamma_inv_trig);
	gmp_printf("Output cross-track error: %Ff\n", debug_1f);
	// END DEBUG

	// Update integral?
	mpz_powm(tmp_x, c_yep, delta_t, N);
	mpz_mul(c_integral, c_integral, tmp_x);
	mpz_mod(c_integral, c_integral, N);

	// Compute c_psi_d
	mpz_powm(tmp_x, c_yep, kp_bar, N);
	mpz_invert(tmp_x, tmp_x, N);
	mpz_powm(inv, c_integral, ki_bar, N);
	mpz_invert(inv, inv, N);
	mpz_mul(c_psi_d, current_c_pi, tmp_x); // Get ciphertext c_pi!
	mpz_mod(c_psi_d, c_psi_d, N);
	mpz_mul(c_psi_d, c_psi_d, inv);
	mpz_mod(c_psi_d, c_psi_d, N);

	// DECRYPT AND DEBUG //
	joye_libert_decrypt(debug_1, c_xn, p, y, msgsize);
	
	rho_inv(debug_1f, debug_1, gamma);
	gmp_printf("c_xn: %Ff\n", debug_1f);

	joye_libert_decrypt(debug_1, c_yn, p, y, msgsize);
	rho_inv(debug_1f, debug_1, gamma);
	gmp_printf("c_yn: %Ff\n", debug_1f);
	//END DEBUG

	// Compute encrypted along-track and cross-track distance to the next waypoint
	// First part
	mpz_t x_tilde;
	mpz_init(x_tilde);
	mpz_invert(inv, c_xn, N);
	//mpz_invert(inv, next_wp_x, N);
	mpz_mul(x_tilde, next_wp_x, inv);
	mpz_mod(x_tilde, x_tilde, N);
	mpz_powm(c_yep, x_tilde, current_sin_pi, N);
	mpz_powm(c_xep, x_tilde, current_cos_pi, N);
	mpz_invert(c_yep, c_yep, N);

	// DECRYPT AND DEBUG //
	joye_libert_decrypt(debug_1, x_tilde, p, y, msgsize);
	
	//std::cout << "Iterate: After debug decrypt\n";
	
	rho_inv(debug_1f, debug_1, gamma);
	gmp_printf("x_tilde: %Ff\n", debug_1f);
	gmp_printf("current_cos_pi: %Zd\n", current_cos_pi);
	gmp_printf("current_sin_pi: %Zd\n", current_sin_pi);
	// //END DEBUG
	
	// Second part
	mpz_invert(inv, next_wp_y, N);
	mpz_mul(tmp_y, c_yn, inv);
	mpz_mod(tmp_y, tmp_y, N);
	mpz_powm(tmp_x, tmp_y, current_sin_pi, N);
	mpz_powm(tmp_y, tmp_y, current_cos_pi, N);

	// Homomorphically add first and second part
	mpz_mul(c_xe, c_xep, tmp_x);
	mpz_mul(c_ye, c_yep, tmp_y);
	mpz_mod(c_xe, c_xe, N);
	mpz_mod(c_ye, c_ye, N);

	// PRINT TO DEBUG HERE
	
	// DECRYPT AND DEBUG //
	joye_libert_decrypt(debug_1, c_xe, p, y, msgsize);
	
	//std::cout << "Iterate: After debug decrypt\n";
	
	rho_inv(debug_1f, debug_1, gamma_inv_trig);
	gmp_printf("c_xe: %Ff\n", debug_1f);
	
	joye_libert_decrypt(debug_1, c_ye, p, y, msgsize);
	rho_inv(debug_1f, debug_1, gamma_inv_trig);
	gmp_printf("c_ye: %Ff\n", debug_1f);
	// //END DEBUG

	/*mpz_invert(inv, c_xp, N);
      	mpz_mul(c_xe, next_wp_x, inv);
	mpz_mod(c_xe, c_xe, N);
	mpz_invert(inv, c_yp, N);
	mpz_mul(c_ye, next_wp_y, inv);
	mpz_mod(c_ye, c_ye, N);
	*/
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

	if (delta < threshold)
	{
		gmp_printf("x_e_bar: %Zd\n", x_e_bar);
		gmp_printf("y_e_bar: %Zd\n", y_e_bar);
		gmp_printf("x_e_p: %Ff\n", x_e_p);
		gmp_printf("y_e_p: %Ff\n", y_e_p);
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
	//std::cout << "Inside rho_inv\n";
	
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
	//std::cout << "Rho_inv: After logic\n";

	mpf_t gamma_f;
	mpf_init(gamma_f);

	mpf_set_z(gamma_f, gamma);
	mpf_set_z(out, in);
	mpf_div(out, out, gamma_f);
	//std::cout << "Rho_inv: Finished\n";
}

