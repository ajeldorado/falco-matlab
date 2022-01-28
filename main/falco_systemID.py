import numpy as np
import tensorflow as tf
import matplotlib.pyplot as plt
import scipy.io as sio

plt.ion()

class SSM(object):
	# linear state space model for the optical system
	# the incoherent light is not considered here, however, it should be easy to include them
	def __init__(self, Jacobian1, Jacobian2, Q0, Q1, R0, R1, R2, n_observ=4):
		# Jacobian matrices
		self.G1_real = tf.Variable(Jacobian1.real, dtype=tf.float64)
		self.G1_imag = tf.Variable(Jacobian1.imag, dtype=tf.float64)
		self.G2_real = tf.Variable(Jacobian2.real, dtype=tf.float64)
		self.G2_imag = tf.Variable(Jacobian2.imag, dtype=tf.float64)
		
		# noise parameters, all of them should be positive
		self.q0 = tf.Variable(np.log(Q0), dtype=tf.float64)
		self.q1 = tf.Variable(np.log(Q1), dtype=tf.float64)
		self.r0 = tf.Variable(np.log(R0), dtype=tf.float64)
		self.r1 = tf.Variable(np.log(R1), dtype=tf.float64)
		self.r2 = tf.Variable(np.log(R2), dtype=tf.float64)

		self.Q0 = tf.exp(self.q0)
		self.Q1 = tf.exp(self.q1)
		self.R0 = tf.exp(self.r0)
		self.R1 = tf.exp(self.r1)
		self.R2 = tf.exp(self.r2)
		self.num_pix = Jacobian1.shape[0]
		self.num_act = Jacobian1.shape[1]
		self.n_observ = int(n_observ)
	def transition(self, Enp, u1, u2):
		# state transition model
		Enp_next = Enp + tf.cast(tf.tensordot(u1, self.G1_real, axes=[[-1], [1]]) + tf.tensordot(u2, self.G2_real, axes=[[-1], [1]]), tf.complex128) + \
					+ 1j * tf.cast(tf.tensordot(u1, self.G1_imag, axes=[[-1], [1]]) + tf.tensordot(u2, self.G2_imag, axes=[[-1], [1]]), tf.complex128)
		return Enp_next
	def transition_covariance(self, Enp, u1, u2):
		# covariance of process/transition noises
		u1_square = tf.reduce_sum(tf.abs(u1)**2, axis=1)
		u2_square = tf.reduce_sum(tf.abs(u2)**2, axis=1)
		Qco = tf.tensordot(tf.expand_dims(u1_square, 1), tf.expand_dims(self.Q1*tf.ones(self.num_pix, dtype=tf.float64), 0), axes=[[1], [0]]) + \
			tf.tensordot(tf.expand_dims(u2_square, 1), tf.expand_dims(self.Q1*tf.ones(self.num_pix, dtype=tf.float64), 0), axes=[[1], [0]]) + self.Q0 + 1e-14
		return Qco
	def observation(self, Enp, u_p1, u_p2):
		# observation model
		n_probe = self.n_observ
		E_p_list = []
		for k in range(n_probe):
			E_p = self.transition(Enp, u_p1[:, k, :], u_p2[:, k, :])
			E_p_list.append(tf.expand_dims(E_p, 1))
		E_p_obs = tf.concat(E_p_list, axis=1)
		I_p_obs = tf.abs(E_p_obs)**2
		return E_p_obs, I_p_obs
	def observation_covariance(self, Enp, u_p1, u_p2):
		# covariance of observation noises
		_, I_p_obs = self.observation(Enp, u_p1, u_p2)
		contrast_p = tf.reduce_mean(I_p_obs, axis=2)

		R = tf.tensordot(tf.expand_dims(self.R0 + self.R1*contrast_p + self.R2*contrast_p**2, axis=-1), tf.ones((1, self.num_pix), dtype=tf.float64), axes=[[-1], [0]]) + 1e-24	
		return R
	def get_params(self):
		# get the parameters of the identified system
		return [self.G1_real, self.G1_imag, self.G2_real, self.G2_imag,
				self.q0, self.q1, self.r0, self.r1, self.r2]


def LSEnet(model, Ip, u1p, u2p):
	# computation graph that defines least squared estimation of the electric field
	delta_Ep_pred = tf.cast(tf.tensordot(u1p, model.G1_real, axes=[[-1], [1]]) + tf.tensordot(u2p, model.G2_real, axes=[[-1], [1]]), tf.complex128) + \
			+ 1j * tf.cast(tf.tensordot(u1p, model.G1_imag, axes=[[-1], [1]]) + tf.tensordot(u2p, model.G2_imag, axes=[[-1], [1]]), tf.complex128)
	delta_Ep_expand = tf.expand_dims(delta_Ep_pred, 2)
	delta_Ep_expand_diff = delta_Ep_expand[:, 1::2, :, :] - delta_Ep_expand[:, 2::2, :, :]
	y = tf.transpose(Ip[:, 1::2, :]-Ip[:, 2::2, :], [0, 2, 1])
	H = tf.concat([2*tf.real(delta_Ep_expand_diff), 2*tf.imag(delta_Ep_expand_diff)], axis=2)
	H = tf.transpose(H, [0, 3, 1, 2])
	Ht_H = tf.matmul(tf.transpose(H, [0, 1, 3, 2]), H)
	Ht_H_inv_Ht = tf.matmul(tf.matrix_inverse(Ht_H+tf.eye(2, dtype=tf.float64)*1e-12), tf.transpose(H, [0, 1, 3, 2]))
	x_new = tf.squeeze(tf.matmul(Ht_H_inv_Ht, tf.expand_dims(y, -1)), -1)
	
	n_observ = model.n_observ
	contrast_p = tf.reduce_mean(Ip, axis=2)

	Rp = tf.tensordot(tf.expand_dims(model.R0 + model.R1*contrast_p + model.R2*contrast_p**2, axis=-1), 
						tf.ones((1, model.num_pix), dtype=tf.float64), axes=[[-1], [0]]) + 1e-24	
	Rp = tf.transpose(Rp, [0, 2, 1])
	R_diff = Rp[:, :, 1::2]+Rp[:, :, 2::2]
	R = tf.matrix_set_diag(tf.concat([tf.expand_dims(tf.zeros_like(R_diff), -1)]*(n_observ//2), -1), R_diff)
	P_new = tf.matmul(tf.matmul(Ht_H_inv_Ht, R), tf.transpose(Ht_H_inv_Ht, [0, 1, 3, 2]))
	Enp_pred_new = tf.cast(x_new[:, :, 0], dtype=tf.complex128) + 1j * tf.cast(x_new[:, :, 1], dtype=tf.complex128)
	return Enp_pred_new, P_new

def KFnet(model, Ip, Enp_old, P_old, u1c, u2c, u1p, u2p):
	# computation graph that defines Kalman filtering estimation of the electric field
	delta_Ep_pred = tf.cast(tf.tensordot(u1p, model.G1_real, axes=[[-1], [1]]) + tf.tensordot(u2p, model.G2_real, axes=[[-1], [1]]), tf.complex128) + \
			+ 1j * tf.cast(tf.tensordot(u1p, model.G1_imag, axes=[[-1], [1]]) + tf.tensordot(u2p, model.G2_imag, axes=[[-1], [1]]), tf.complex128)
	delta_Ep_expand = tf.expand_dims(delta_Ep_pred, 2)
	delta_Ep_expand_diff = delta_Ep_expand[:, 1::2, :, :] - delta_Ep_expand[:, 2::2, :, :]
	y = tf.transpose(Ip[:, 1::2, :]-Ip[:, 2::2, :], [0, 2, 1])
	H = tf.concat([2*tf.real(delta_Ep_expand_diff), 2*tf.imag(delta_Ep_expand_diff)], axis=2)
	H = tf.transpose(H, [0, 3, 1, 2])

	n_observ = model.n_observ
	contrast_p = tf.reduce_mean(Ip, axis=2)
	Rp = tf.tensordot(tf.expand_dims(model.R0 + model.R1*contrast_p + model.R2*contrast_p**2, axis=-1), 
						tf.ones((1, model.num_pix), dtype=tf.float64), axes=[[-1], [0]]) + 1e-24	
	Rp = tf.transpose(Rp, [0, 2, 1])
	R_diff = Rp[:, :, 1::2]+Rp[:, :, 2::2]
	R = tf.matrix_set_diag(tf.concat([tf.expand_dims(tf.zeros_like(R_diff), -1)]*(n_observ//2), -1), R_diff)

	Qco = model.transition_covariance(Enp_old, u1c, u2c)
	Q = tf.concat([tf.expand_dims(Qco, 2), tf.expand_dims(Qco, 2)], axis=2)
	Q = tf.matrix_set_diag(tf.concat([tf.expand_dims(tf.zeros_like(Q), -1)]*2, -1), Q)

	# state and covariance prediction
	Enp_pred = model.transition(Enp_old, u1c, u2c)
	Enp_expand = tf.expand_dims(Enp_pred, 2)
	x = tf.concat([tf.real(Enp_expand), tf.imag(Enp_expand)], axis=2)
	P = P_old + Q

	# state and covariance update
	y_model = tf.squeeze(tf.matmul(H, tf.expand_dims(x, -1)), -1)
	S = R + tf.matmul(tf.matmul(H, P), tf.transpose(H, [0, 1, 3, 2]))
	K = tf.matmul(tf.matmul(P, tf.transpose(H, [0, 1, 3, 2])), tf.matrix_inverse(S))
	x_new = x + tf.squeeze(tf.matmul(K, tf.expand_dims(y-y_model, -1)), -1)
	P_new = P - tf.matmul(tf.matmul(K, H), P)
	
	Enp_pred_new = tf.cast(x_new[:, :, 0], dtype=tf.complex128) + 1j * tf.cast(x_new[:, :, 1], dtype=tf.complex128)
	return Enp_pred_new, P_new



# if __name__ == "__main__":
def linear_vl(Q0=1e-10, Q1=1e-7, R0=1e-14, R1=1e-9, R2=1e-2, 
				lr=1e-7, lr2=1e-2, epoch=10, print_flag=False, path2data='/Users/ajriggs/Repos/falco-matlab/data/jac/'):
	# load training data and biased Jacobian matrices
	data_train = sio.loadmat( (path2data+'data_train.mat'))
	u1_train = data_train['data_train']['u1'][0, 0]
	u2_train = data_train['data_train']['u2'][0, 0]
	u1p_train = data_train['data_train']['u1p'][0, 0]
	u2p_train = data_train['data_train']['u2p'][0, 0]
	image_train = data_train['data_train']['I'][0, 0]

	Jacobian = sio.loadmat((path2data+'jacStruct.mat'))
	G1 = Jacobian['jacStruct']['G1'][0, 0]
	G2 = Jacobian['jacStruct']['G2'][0, 0]

	# extract the parameters of the system
	n_step = u1_train.shape[1] # number of control steps
	n_act = u1p_train.shape[0] # number of active actuators on the DM
	n_pix = image_train.shape[0] # number of pixels in the dark hole
	n_pair = u1p_train.shape[1]/2 # number of probing pairs in each control step
	n_image = 2 * n_pair + 1 # number of probe images in each control step

	u1p_train = np.concatenate([np.zeros((n_act, 1, n_step)), u1p_train], axis=1)
	u2p_train = np.concatenate([np.zeros((n_act, 1, n_step)), u2p_train], axis=1)

	# define the placeholders for the computation graph
	mc_sampling = 5
	u1c = tf.placeholder(tf.float64, shape=(None, n_act))
	u2c = tf.placeholder(tf.float64, shape=(None, n_act))
	u1p = tf.placeholder(tf.float64, shape=(None, n_image, n_act))
	u2p = tf.placeholder(tf.float64, shape=(None, n_image, n_act))
	Enp_old = tf.placeholder(tf.complex128, shape=(None, n_pix))
	Iinco_old = tf.placeholder(tf.float64, shape=(None, n_pix))
	Iinco = tf.placeholder(tf.float64, shape=(None, n_pix))
	Ip = tf.placeholder(tf.float64, shape=(None, n_image, n_pix))
	P_old = tf.placeholder(tf.float64, shape=(None, n_pix, 2, 2))
	noises = tf.placeholder(tf.float64, shape=(None, n_pix, 2, mc_sampling))
	learning_rate = tf.placeholder(tf.float64, shape=())
	learning_rate2 = tf.placeholder(tf.float64, shape=())


	# define the optical model as a state space model (SSM) or a neural network
	# parameters Q0, Q1, R0, R1, R2 define the noise covariance of the state space model
	# process noises covariance: Q = Q0 + Q1 * sum(uc^2)
	# observation noises covariance: R = R0 + R1 * probe_contrast + R2 * probe_contrast^2
	model = SSM(G1, G2, Q0, Q1, R0, R1, R2, n_image)

	# define the relations of the control/probe inputs, camera images and hidden electric fields
	Enp_pred = model.transition(Enp_old, u1c, u2c)
	Qco = model.transition_covariance(Enp_old, u1c, u2c)

	# Enp_est, P_est = LSEnet(model, Ip, u1p, u2p)
	Enp_est2, P_est2 = LSEnet(model, Ip, u1p, u2p)
	Enp_est, P_est = KFnet(model, Ip, Enp_old, P_old, u1c, u2c, u1p, u2p)

	P_est_sqrt = tf.cholesky(P_est+tf.eye(2, dtype=tf.float64)*1e-14)
	Enp_est_noises = tf.matmul(P_est_sqrt, noises)
	Enp_est_noises_complex = tf.cast(Enp_est_noises[:, :, 0, :], tf.complex128) + 1j * tf.cast(Enp_est_noises[:, :, 1, :], tf.complex128)

	Enp_est_samples = Enp_est + Enp_est_noises_complex[:, :, 0]
	_, Ip_pred = model.observation(Enp_est_samples, u1p, u2p)
	Rp = model.observation_covariance(Enp_est_samples, u1p, u2p)
	for k in range(1, mc_sampling):
		Enp_est_samples = Enp_est + Enp_est_noises_complex[:, :, k]
		_, Ip_pred_now = model.observation(Enp_est_samples, u1p, u2p)
		Ip_pred += Ip_pred_now
		Rp += model.observation_covariance(Enp_est_samples, u1p, u2p)
	Ip_pred /= mc_sampling
	Rp /= mc_sampling

	# evidence lower bound (elbo): cost function for system identification
	# we need to maximize the elbo for system ID
	elbo = - tf.reduce_sum(tf.abs(Ip-Ip_pred)**2 / Rp) - tf.reduce_sum(tf.log(2*np.pi*Rp)) - \
			(tf.reduce_sum(tf.abs(Enp_pred-Enp_est)**2 / Qco) + tf.reduce_sum(2 * tf.log(Qco)) - \
			 tf.reduce_sum(tf.linalg.logdet(P_est)) + tf.reduce_sum(tf.trace(P_est) / Qco))

	# mean squared error (MSE): a metric for checking the system ID results
	MSE = tf.reduce_sum(tf.abs(Ip - Ip_pred)**2)

	params_list = model.get_params() # parameters to be identified


	# start identifying/learning the model parameters
	train_Jacobian = tf.train.AdamOptimizer(learning_rate=learning_rate, 
											beta1=0.99, beta2=0.9999, epsilon=1e-08).minimize(-elbo, var_list=params_list[0:4])
	train_noise_coef = tf.train.AdamOptimizer(learning_rate=learning_rate2, 
											beta1=0.99, beta2=0.9999, epsilon=1e-08).minimize(-elbo, var_list=params_list[4::])
	train_op = tf.group(train_Jacobian, train_noise_coef)
	# train_op = tf.train.AdamOptimizer(learning_rate=learning_rate, 
	# 										beta1=0.99, beta2=0.9999, epsilon=1e-08).minimize(-elbo, var_list=params_list)
	init = tf.global_variables_initializer()
	elbo_list = []
	mse_list = []
	epoch = int(epoch)
	with tf.Session() as sess:
		sess.run(init)
		# mse = sess.run(MSE, feed_dict={Ip: np.transpose(image_train, [2, 1, 0]),
		# 								u1p: np.transpose(u1p_train, [2, 1, 0]),
		# 								u2p: np.transpose(u2p_train, [2, 1, 0]),
		# 								noises: np.zeros((n_step, n_pix, 2,	mc_sampling))})

		Enp_est_values, P_est_values = sess.run([Enp_est2, P_est2], feed_dict={Ip: np.transpose(image_train, [2, 1, 0]),
														u1p: np.transpose(u1p_train, [2, 1, 0]),
														u2p: np.transpose(u2p_train, [2, 1, 0])})
		mse = sess.run(MSE, feed_dict={Ip: np.transpose(image_train[:, :, 1:n_step], [2, 1, 0]),
										u1p: np.transpose(u1p_train[:, :, 1:n_step], [2, 1, 0]),
										u2p: np.transpose(u2p_train[:, :, 1:n_step], [2, 1, 0]),
										u1c: np.transpose(u1_train[:,0:n_step-1]), 
										u2c: np.transpose(u2_train[:,0:n_step-1]),
										Enp_old: Enp_est_values[0:n_step-1, :],
										P_old: P_est_values[0:n_step-1, :, :, :],
										noises: np.zeros((n_step-1, n_pix, 2, mc_sampling))})

		mse_list.append(mse)
		print('initial MSE: {}'.format(mse))

		for k in range(epoch):
			# Enp_est_values = sess.run(Enp_est, feed_dict={Ip: np.transpose(image_train, [2, 1, 0]),
			# 											u1p: np.transpose(u1p_train, [2, 1, 0]),
			# 											u2p: np.transpose(u2p_train, [2, 1, 0])})

			# sess.run(train_op, feed_dict={noises: np.random.normal(size=(n_step-1, n_pix, 2, mc_sampling)),
			# 								Enp_old: Enp_est_values[0:n_step-1, :],
			# 								Ip: np.transpose(image_train[:, :, 1:n_step], [2, 1, 0]),
			# 								u1c: np.transpose(u1_train[:,0:n_step-1]), u2c: np.transpose(u2_train[:,0:n_step-1]),
			# 								u1p: np.transpose(u1p_train[:, :, 1:n_step], [2, 1, 0]), u2p: np.transpose(u2p_train[:, :, 1:n_step], [2, 1, 0]),
			# 								learning_rate: lr, learning_rate2: lr2})
			# mse = sess.run(MSE, feed_dict={Ip: np.transpose(image_train, [2, 1, 0]),
			# 							u1p: np.transpose(u1p_train, [2, 1, 0]),
			# 							u2p: np.transpose(u2p_train, [2, 1, 0]),
			# 							noises: np.zeros((n_step, n_pix, 2, mc_sampling))})

			Enp_est_values0, P_est_values0 = sess.run([Enp_est2, P_est2], feed_dict={Ip: np.transpose(np.expand_dims(image_train[:, :, 0], -1), [2, 1, 0]),
														u1p: np.transpose(np.expand_dims(u1p_train[:, :, 0], -1), [2, 1, 0]),
														u2p: np.transpose(np.expand_dims(u2p_train[:, :, 0], -1), [2, 1, 0])})

			Enp_est_values[0, :] = np.squeeze(Enp_est_values0)
			P_est_values[0, :, :, :] = np.squeeze(P_est_values0)
			for i in range(1, n_step):
				Enp_est_values_now, P_est_values_now = sess.run([Enp_est, P_est], feed_dict={Ip: np.transpose(np.expand_dims(image_train[:, :, i], -1), [2, 1, 0]),
											u1p: np.transpose(np.expand_dims(u1p_train[:, :, i], -1), [2, 1, 0]),
											u2p: np.transpose(np.expand_dims(u2p_train[:, :, i], -1), [2, 1, 0]),
											u1c: np.transpose(np.expand_dims(u1_train[:,i-1], -1)), 
											u2c: np.transpose(np.expand_dims(u2_train[:,i-1], -1)),
											Enp_old: np.expand_dims(Enp_est_values[i-1, :], 0),
											P_old: np.expand_dims(P_est_values[i-1, :, :, :], 0),
											noises: np.zeros((1, n_pix, 2, mc_sampling))})
				Enp_est_values[i, :] = np.squeeze(Enp_est_values_now)
				P_est_values[i, :, :, :] = np.squeeze(P_est_values_now)

			sess.run(train_op, feed_dict={noises: np.random.normal(size=(n_step-1, n_pix, 2, mc_sampling)),
											Enp_old: Enp_est_values[0:n_step-1, :], P_old: P_est_values[0:n_step-1, :, :, :],
											Ip: np.transpose(image_train[:, :, 1:n_step], [2, 1, 0]),
											u1c: np.transpose(u1_train[:,0:n_step-1]), u2c: np.transpose(u2_train[:,0:n_step-1]),
											u1p: np.transpose(u1p_train[:, :, 1:n_step], [2, 1, 0]), u2p: np.transpose(u2p_train[:, :, 1:n_step], [2, 1, 0]),
											learning_rate: lr, learning_rate2: lr2})


			Enp_est_values, P_est_values = sess.run([Enp_est2, P_est2], feed_dict={Ip: np.transpose(image_train, [2, 1, 0]),
														u1p: np.transpose(u1p_train, [2, 1, 0]),
														u2p: np.transpose(u2p_train, [2, 1, 0])})
			mse = sess.run(MSE, feed_dict={Ip: np.transpose(image_train[:, :, 1:n_step], [2, 1, 0]),
										u1p: np.transpose(u1p_train[:, :, 1:n_step], [2, 1, 0]),
										u2p: np.transpose(u2p_train[:, :, 1:n_step], [2, 1, 0]),
										u1c: np.transpose(u1_train[:,0:n_step-1]), 
										u2c: np.transpose(u2_train[:,0:n_step-1]),
										Enp_old: Enp_est_values[0:n_step-1, :],
										P_old: P_est_values[0:n_step-1, :, :, :],
										noises: np.zeros((n_step-1, n_pix, 2, mc_sampling))})


			if mse < mse_list[-1]:
				mse_list.append(mse)
			else:
				break
			if print_flag:
				print('epoch {} MSE: {}'.format(k, mse))
				# print('Q0: {}, Q1: {}, R0: {}, R1: {}, R2: {}'.format(sess.run(model.Q0), sess.run(model.Q1),
				# 													sess.run(model.R0), sess.run(model.R1), sess.run(model.R2)))

		params_values_list = []
		for param in params_list:
			params_values_list.append(sess.run(param))

	sio.savemat(path2data+'jacStructLearned.mat', {'G1': params_values_list[0]+1j*params_values_list[1],
										'G2': params_values_list[2]+1j*params_values_list[3],
										'noise_coef': np.array(params_values_list[-5::])})