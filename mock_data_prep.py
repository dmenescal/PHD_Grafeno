import random 
import numpy as np
from numpy import genfromtxt


def load_train_test_mockset(source_file='char8_12x10.txt', split_ratio=0.8):
    
    def img_noise(sample, prob):
        noisy_img = np.zeros(sample.shape)
        noisy_samples = np.zeros(sample.shape)
        threshold = 1 - prob
        for i in range(sample.shape[0]):
            for j in range(sample.shape[1]):
                rdn = random.random()
                if (rdn < prob):
                    noisy_img[i][j] = 0
                elif (rdn > threshold):
                    noisy_img[i][j] = 1
                    noisy_samples[i][j] = 1
                else:
                    noisy_img[i][j] = sample[i][j]
        return noisy_img
    
    data = genfromtxt(source_file, delimiter=',')

    char_num = data.shape[1]

    char_list = []

    for i in range(char_num):
        char_vec_img = data[:, i]
        char_mat_img = char_vec_img.reshape((12, 10))
        char_list.append(char_mat_img)

    char_data_set = np.array(char_list)

    norm_dataset = (char_data_set -np.min(char_data_set))/(np.max(char_data_set)-np.min(char_data_set))

    # Geração de um conjunto de probabilidades para aplicar o ruido na imagem
    noise_prob_levels = np.arange(0, 0.10, 0.0005)

    # Construção de uma base de amostras ruidosas
    label_samples = []
    for label in range(8):
        samples = []
        for sample in range(len(noise_prob_levels)):
            noisy_img = img_noise(norm_dataset[label], noise_prob_levels[i])
            samples.append(noisy_img)
        label_samples.append(samples)

    sample_set = np.array(label_samples)
    
    # Separação da base de amostras em base de treino e base de teste, na proporção 80-20
    train_lower_idx = 0
    train_upper_idx = int(sample_set.shape[1]*split_ratio)
    test_lower_idx = int(sample_set.shape[1]*split_ratio)
    test_upper_idx = test_lower_idx + int(sample_set.shape[1]*(1-split_ratio)) + 1

    train_set = sample_set[:, train_lower_idx:train_upper_idx]
    test_set = sample_set[:, test_lower_idx:test_upper_idx]
    
    return train_set, test_set
