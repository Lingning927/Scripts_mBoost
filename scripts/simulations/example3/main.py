import torch
import torch.nn as nn
import pandas as pd
from torch.utils.data import TensorDataset, DataLoader

# 数据加载和预处理
def load_data(x_path, y_path):
    X = pd.read_csv(x_path).values
    Y = pd.read_csv(y_path).values
    X_tensor = torch.tensor(X, dtype=torch.float32)
    Y_tensor = torch.tensor(Y, dtype=torch.float32)
    dataset = TensorDataset(X_tensor, Y_tensor)
    return DataLoader(dataset, batch_size=100, shuffle=True)

# 神经网络模型
class NeuralNetwork(nn.Module):
    def __init__(self, input_size, output_size, num_layers, layer_size):
        super(NeuralNetwork, self).__init__()
        layers = [nn.Linear(input_size, layer_size), nn.ReLU()]
        for _ in range(num_layers - 1):
            layers.append(nn.Linear(layer_size, layer_size))
            layers.append(nn.ReLU())
        layers.append(nn.Linear(layer_size, output_size))
        self.model = nn.Sequential(*layers)

    def forward(self, x):
        return self.model(x)

# 训练模型
def train_model(model, data_loader, epochs=100):
    criterion = nn.MSELoss()
    optimizer = torch.optim.Adam(model.parameters())
    for epoch in range(epochs):
        for inputs, targets in data_loader:
            optimizer.zero_grad()
            outputs = model(inputs)
            loss = criterion(outputs, targets)
            loss.backward()
            optimizer.step()
 
# 主函数
def main():
    x_path = 'scripts/simulations/example3/X_train.csv'
    y_path = 'scripts/simulations/example3/Y_train.csv'
    x_test_path = 'scripts/simulations/example3/X_test.csv'
    data_loader = load_data(x_path, y_path)

    input_size = 784  # 根据您的数据集调整
    output_size = 1   # 根据您的数据集调整

    depths = [2, 3, 4, 5, 6, 7]
    widths = [2, 4, 8, 16, 32, 64, 128]
    predictions_matrix_train = []
    predictions_matrix_test = []

    X_train = pd.read_csv(x_path).values
    X_train_tensor = torch.tensor(X_train, dtype=torch.float32)

    X_test = pd.read_csv(x_test_path).values
    X_test_tensor = torch.tensor(X_test, dtype=torch.float32)

    for depth in depths:
        for width in widths:
            model = NeuralNetwork(input_size, output_size, depth, width)
            train_model(model, data_loader)

            with torch.no_grad():
                predictions_train = model(X_train_tensor).numpy().flatten()
                predictions_test = model(X_test_tensor).numpy().flatten()
            predictions_matrix_train.append(predictions_train)
            predictions_matrix_test.append(predictions_test)

    predictions_df_train = pd.DataFrame(predictions_matrix_train).T
    predictions_df_train.columns = [f'Depth {d}_Width {w}' for d in depths for w in widths]
    predictions_df_train.to_csv('scripts/simulations/example3/predictions_matrix_train.csv', index=False)

    predictions_df_test = pd.DataFrame(predictions_matrix_test).T
    predictions_df_test.columns = [f'Depth {d}_Width {w}' for d in depths for w in widths]
    predictions_df_test.to_csv('scripts/simulations/example3/predictions_matrix_test.csv', index=False)

if __name__ == '__main__':
    main()
