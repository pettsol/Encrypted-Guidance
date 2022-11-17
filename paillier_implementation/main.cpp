#include <iostream>
#include <fstream>
#include <unistd.h>
#include <chrono>

#include "paillier/paillier.h"
#include "he_guidance.h"
#include "encoder.h"


// HELPER FUNCTION FOR LOGGING PURPOSES
void rho_inv(mpf_t out, mpz_t in, mpz_t gamma)
{
	uint32_t msgsize = 32;	
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


int main()
{
	size_t ct_size_byte = 256;
	uint32_t keysize = 2048;
	uint32_t msgsize = 32;

	gmp_randstate_t rand_state;
	gmp_randinit_mt(rand_state);

	mpz_t N, N2, lambda, g, p2, q2, mu, gamma_p,
	      gamma_kp, gamma_ki, gamma_inverse, gamma_inv_trig, gamma_time;
	mpz_init(N);
	mpz_init(N2);
	mpz_init(lambda);
	mpz_init(g);
	mpz_init(p2);
	mpz_init(q2);
	mpz_init(mu);
	mpz_init_set_ui(gamma_p, 100);
	mpz_init_set_ui(gamma_kp, 10000);
	mpz_init_set_ui(gamma_ki, 10000);
	mpz_init_set_ui(gamma_inverse, 100000000);
	mpz_init_set_ui(gamma_inv_trig, 10000);
	mpz_init_set_ui(gamma_time, 1);

	paillier_keygen(N, N2, p2, q2, g, lambda, mu, rand_state, keysize);

	float kp = 0.1;
	//float ki = 1;
	float ki = 0.001;
	float delta_t = 1;
	float threshold = 1;

	Encrypted_ilos_guidance state(N, N2, lambda, g, p2, q2, mu, gamma_p, gamma_kp, gamma_ki, gamma_inverse, gamma_inv_trig, gamma_time, threshold, kp, ki, delta_t, msgsize);


	// Create two waypoints (1,1) and (2,2)
	float waypoints[12] = {0, 0, 10, 10, 15, 36, 30, 46, 58, 40};

	// Encrypt the waypoints
	state.preprocessing(waypoints, 5);
	
	// Pass the USV NED position (0,0) to quantize and encrypt
	uint32_t b = 0;
	float x_pos = -4;
	float y_pos = 2;
	float speed = 1;
	mpz_t c_x, c_y;
	mpz_init(c_x);
	mpz_init(c_y);
	
	mpf_t psi_d_f;
	mpf_init(psi_d_f);
	
	mpz_t c_psi_d, c_xe, c_ye;
	mpz_init(c_psi_d);
	mpz_init(c_xe);
	mpz_init(c_ye);

	float heading = 0;
	float heading_rate = 0;
	float desired_heading = 0;
	float kp_yaw = 1;

	float K, T, delta;
	K = 1/41.4; T = 1;

	// Open log - Log both position and yaw?
	std::ofstream log_position("position.txt");
	std::ofstream log_yaw("yaw.txt");
	std::ofstream log_desired_yaw("desired_yaw.txt");
	std::ofstream log_waypoints("waypoints.txt");
	std::ofstream log_cte("cte.txt");

	// Log time
	std::ofstream log_encryption("encryption_delay.txt");
	std::ofstream log_evaluation("evaluation_delay.txt");
	std::ofstream log_decryption("decryption_delay.txt");

	for (int i = 0; i < 5; i++)
	{
		log_waypoints << waypoints[2*i] << "\t" << waypoints[2*i+1] << std::endl;
	}

	uint8_t count = 0;
	while (count < 5)
	{
		// LOG ENCRYPTION TIME
		auto start_time = std::chrono::system_clock::now();
		state.quantize_and_encrypt(c_x, c_y, x_pos, y_pos);
		auto stop_time = std::chrono::system_clock::now();
		auto latency = std::chrono::duration_cast<std::chrono::microseconds>(stop_time - start_time);
		log_encryption << latency.count() << std::endl;
		// END OF ENCRYPTION TIME LOG

		// LOG ENCRYPTED GUIDANCE TIME
		// Pass the encrypted NED position to the controller
		start_time = std::chrono::system_clock::now();
		state.iterate(c_psi_d, c_xe, c_ye, &b, c_x, c_y);
		stop_time = std::chrono::system_clock::now();
		latency = std::chrono::duration_cast<std::chrono::microseconds>(stop_time - start_time);
		log_evaluation << latency.count() << std::endl;

		// END OF ENCRYPTED GUIDANCE TIME LOG
		//
		// LOG DECRYPTION TIME
		// Decrypt and recover the desired heading
		start_time = std::chrono::system_clock::now();
		state.decrypt_and_recover(psi_d_f, &b, c_psi_d, c_xe, c_ye);
		stop_time = std::chrono::system_clock::now();
		latency = std::chrono::duration_cast<std::chrono::microseconds>(stop_time - start_time);
		log_decryption << latency.count() << std::endl;
		// END OF DECRYPTION TIME LOG
		if (b)
		{
			count++;
		}

		// Print out the desired heading
		desired_heading = mpf_get_d(psi_d_f);

		std::cout << "Desired heading: " << desired_heading << std::endl;

		
		// Simulate dynamic system? Use a first order model for heading? What about position?
		// Constant speed, varying heading? Set surge speed to 0.1 m/s, rest to 0? Use 
		// desired heading directly?
		for (int i = 0; i < 5; i++)
		{
			// Proportional heading control with first-order Nomoto model
			// We control yaw-rate!
			heading_rate = heading_rate + (-(heading_rate/T) + kp_yaw*(desired_heading - heading))*0.1;
			heading = heading + heading_rate*0.1;
			std::cout << "Heading_rate: " << heading_rate << std::endl;
			std::cout << "Heading: " << heading << std::endl;

			log_desired_yaw << desired_heading << std::endl;
			log_yaw << heading << std::endl;

			x_pos = x_pos + 0.1*speed*cos(heading);
			y_pos = y_pos + 0.1*speed*sin(heading);

			log_position << x_pos << "\t" << y_pos << std::endl;

			std::cout << "x: " << x_pos << std::endl;
			std::cout << "y: " << y_pos << std::endl;
		}

		sleep(0.1);
	}
}
