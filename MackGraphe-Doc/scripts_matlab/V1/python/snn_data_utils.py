import pandas as pd
import numpy as np
import os
from skimage import io, transform
import torch
from torch.utils.data import Dataset, DataLoader

class WaveGuideDataset(Dataset):
    """Light Waveguide dataset."""

    def __init__(self, csv_file, root_dir, transform=None):
        """
        Args:
            csv_file (string): Path to the csv file with annotations.
            root_dir (string): Directory with all the images.
            transform (callable, optional): Optional transform to be applied
                on a sample.
        """
        self.deciles_frame = pd.read_csv(csv_file)
        self.root_dir = root_dir
        self.transform = transform

    def __len__(self):
        return len(self.deciles_frame)

    def __getitem__(self, idx):
        if torch.is_tensor(idx):
            idx = idx.tolist()

        img_name = os.path.join(self.root_dir,
                                self.deciles_frame.iloc[idx, 0]).replace("\\","/")

        image = io.imread(img_name)
        image = image/255
        deciles = self.deciles_frame.iloc[idx, 2:12]
        deciles = np.array([deciles]).astype('float').reshape(-1, 10)
        sample = {'image': image, 'deciles': deciles}

        if self.transform:
            sample = self.transform(sample)

        return sample

class ToTensor(object):
    """Convert ndarrays in sample to Tensors."""

    def __call__(self, sample):
        image, deciles = sample['image'], sample['deciles']

        # swap color axis because
        # numpy image: H x W x C
        # torch image: C x H x W
        image = image.transpose((2, 0, 1))
        return {'image': torch.as_tensor(image),
                'deciles': torch.as_tensor(deciles)}