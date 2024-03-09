import torch
import torch.nn as nn
from collections import OrderedDict
import math
import torch.nn.functional as F
import torch.utils.data
import torch.nn.modules.activation as activation


class ATAC_CNN_Transformer(nn.Module):
    def __init__(self, seq_len, cnn_kernel_size, cnn_nfilter, cnn_maxpooling_size, cnn_dropout,
                 tf_features, tf_nhead, tf_dropout, tf_nlayer, dense_dropout):
        super(ATAC_CNN_Transformer, self).__init__()
        # Define one layer conv1d
        self.conv = nn.Sequential(
            nn.Conv1d(in_channels=4, out_channels=cnn_nfilter, kernel_size=cnn_kernel_size, stride=1),
            nn.BatchNorm1d(cnn_nfilter),
            nn.MaxPool1d(kernel_size=cnn_maxpooling_size, stride=cnn_maxpooling_size),
            nn.Dropout(cnn_dropout),
            nn.ReLU()
        )
        # Define transformer encoding layer
        encoder_layer = nn.TransformerEncoderLayer(d_model=tf_features, nhead=tf_nhead)
        self.encoder = nn.TransformerEncoder(encoder_layer, num_layers=tf_nlayer) 
        # Define dense layer
        self.dense = nn.Sequential(
            # nn.Flatten(),
            nn.Dropout(tf_dropout),
            nn.Linear(tf_features*cnn_nfilter, 200),
            nn.ReLU(),
            nn.Dropout(dense_dropout),
            nn.Linear(200, 24),
            nn.ReLU(),
            nn.Linear(24, 1),
            nn.Sigmoid(),
        )
        
    def forward(self, inputs):
        inputs = inputs.permute(0,2,1)
        cnn_output = self.conv(inputs)
        tf_output = self.encoder(cnn_output)
        tf_output = tf_output.view(tf_output.size(0), -1)
        dense_output = self.dense(tf_output)
        
        return dense_output