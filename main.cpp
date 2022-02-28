#include <iostream>
#include <fstream>
#include <unistd.h>

#include "he_guidance.h"
#include "encoder.h"

int main()
{
	size_t ct_size_byte = 256;
	uint32_t keysize = 2048;
	uint32_t msgsize = 32;

	// Set N, y, p?
	std::string N_string = "C09858A300D8852F23BC728E0585FB5BD8DB7490394F62122B75C1769A642B63136E87DF442CBFF03F12055DC173829247C492517822924D6C5C2390CF71D00282081B23327AA3A1F9C103F5395E3F243FA4D14E3ACA47D530D5D4341A1713CF0425C3F755ED6ED3E0A486DE977DEA6D0A149BBFF66BFCF3BAFB4217B8975A9EA695CB9066FBB00B7C97656DC171E6095898343D5DC98B24D115CA12FF17392DABE07164B3E8CB26890EBCB528E9ADF32597DBFB048A7D6D2494345FE0A4E107D1870395CDEE2B438BAD17C150EA44D4CA8D2290CF4F3113263D8CAF5E6B13366DA2F584D18519FAA7E847973B5F31BA7B970FF831C89201B1AFB44400000001";
        std::string y_string = "53E61608834634EAE18E60C5B991B9F8D71B2D971CBE5AC9E09F4814ADDAB421EFDCC2870D2C92C87003FCFF55CCBA1D4F22F5AB90950FB020F8BE80BA9B4C7CA011F74C2D41581F0036D233B5E8E58B6DD5CA6DB0625D764B927A43FE78844090C6843F29A331B76F8ECE93E7E313ECCB9BCB6ED2330923899AAE43A0FD2430CB6772793755E74862E61E2AC376CFAB9D61827E646421B28E9E0E2ACA4625731AEBBB69EA37E0FA859E499B8A186C8EE6196954170EB8068593F0D764150A6D2E5D3FEA7D9D0D33AC553EECD5C3F27A310115D283E49377820195C8E67781B6F112A625B14B747FA4CC13D06EBA0917246C775F5C732865701AE9349EA8729C";	
	std::string p_string = "DE0BBADE38204E63359A46E672A8D0A2FD5300692AB48F9EF732F5C3FA212B90C98229BBB79BECE734A622154C904DCE9A0F53D4A88B3E558EF7612F6694CE7518F204FE6846AEB6F58174D57A3372363C0D9FCFAA3DC18B1EFF7E89BF7678636580D17DD84A873B14B9C0E1680BBDC87647F3C382902D2F58D6541300000001";	
	
	uint8_t N_array[ct_size_byte] = {0};
        uint8_t y_array[ct_size_byte] = {0};
	uint8_t p_array[ct_size_byte/2] = {0};

        std::cout << "Hex decoding\n";
        hex_decode(N_array, N_string.data(), N_string.size());
        hex_decode(y_array, y_string.data(), y_string.size());
	hex_decode(p_array, p_string.data(), p_string.size());

	mpz_t N, y, p, gamma, gamma_inverse, gamma_inv_trig, gamma_time;
	mpz_init(N);
	mpz_init(y);
	mpz_init(p);
	mpz_init_set_ui(gamma, 10);
	mpz_init_set_ui(gamma_inverse, 1000);
	mpz_init_set_ui(gamma_inv_trig, 100);
	mpz_init_set_ui(gamma_time, 1);

	float kp = 0.25;
	//float ki = 1;
	float ki = 0.05;
	float delta_t = 1;
	float threshold = 1;

        std::cout << "Mpz import\n";
        mpz_import(N, ct_size_byte, 1, 1, 0, 0, N_array);
        mpz_import(y, ct_size_byte, 1, 1, 0, 0, y_array);
	mpz_import(p, ct_size_byte/2, 1, 1, 0, 0, p_array);

	std::cout << "instantiating state\n";
	Encrypted_ilos_guidance state(N, y, p, gamma, gamma_inverse, gamma_inv_trig, gamma_time, threshold, kp, ki, delta_t, msgsize);


	// Create two waypoints (1,1) and (2,2)
	float waypoints[12] = {0, 0, 5, 10, 10, 18, 15, 23, 20, 20};

	// Encrypt the waypoints
	state.preprocessing(waypoints, 5);
	std::cout << "After preprocessing\n";
	
	
	// Pass the USV NED position (0,0) to quantize and encrypt
	uint32_t b = 0;
	float x_pos = -1;
	float y_pos = 0;
	float speed = 0.1;
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
	float kp_yaw = 1.5;

	float K, T, delta;
	K = 0.5; T = 1;

	// Open log - Log both position and yaw?
	std::ofstream log_position("position.txt");
	std::ofstream log_yaw("yaw.txt");
	std::ofstream log_desired_yaw("desired_yaw.txt");
	std::ofstream log_waypoints("waypoints.txt");

	for (int i = 0; i < 5; i++)
	{
		log_waypoints << waypoints[2*i] << "\t" << waypoints[2*i+1] << std::endl;
	}

	uint8_t count = 0;
	while (count < 5)
	{
		std::cout << "outer b: " << b << std::endl;
		state.quantize_and_encrypt(c_x, c_y, x_pos, y_pos, &b);

		//std::cout << "After quantize and encrypt\n";
		// Pass the encrypted NED position to the controller
		state.iterate(c_psi_d, c_xe, c_ye, &b, c_x, c_y);

		//std::cout << "After iterate\n";
		// Decrypt and recover the desired heading
		state.decrypt_and_recover(psi_d_f, &b, c_psi_d, c_xe, c_ye);

		if (b)
		{
			count++;
		}

		//std::cout << "After decrypt and recover\n";
		// Print out the desired heading
		desired_heading = mpf_get_d(psi_d_f);

		std::cout << "Desired heading: " << desired_heading << std::endl;

		
		// Simulate dynamic system? Use a first order model for heading? What about position?
		// Constant speed, varying heading? Set surge speed to 0.1 m/s, rest to 0? Use 
		// desired heading directly?
		for (int i = 0; i < 10; i++)
		{
			// Proportional heading control with first-order Nomoto model
			// We control yaw-rate!
			heading_rate = heading_rate + (-(heading_rate/T) + kp_yaw*(desired_heading - heading))*0.1;
			heading = heading + heading_rate*0.1;
			std::cout << "Heading_rate: " << heading_rate << std::endl;
			std::cout << "Heading: " << heading << std::endl;
			//heading = 0.9*heading + kp_yaw*(desired_heading-heading);

			log_desired_yaw << desired_heading << std::endl;
			log_yaw << heading << std::endl;

			x_pos = x_pos + 0.1*speed*cos(heading);
			y_pos = y_pos + 0.1*speed*sin(heading);

			log_position << x_pos << "\t" << y_pos << std::endl;

			std::cout << "x: " << x_pos << std::endl;
			std::cout << "y: " << y_pos << std::endl;
		}

		//std::cout << "x: " << x_pos << std::endl;
		//std::cout << "y: " << y_pos << std::endl;



		sleep(0.1);
	}
	/*
	// RHO AND RHO INV DEBUG
	std::cout << "TESTING RHO AND RHO INV\n\n";
	mpf_t test_f;
	mpf_init_set_d(test_f, -1.0);
	gmp_printf("test_f: %Ff\n", test_f);

	mpz_t test_z;
	mpz_init(test_z);

	state.rho(test_z, test_f, gamma);
	gmp_printf("test_z: %Zd\n", test_z);

	state.rho_inv(test_f, test_z, gamma);
	gmp_printf("test_f: %Ff\n", test_f);

	// END RHO AND RHO INV DEBUG
	*/
}
