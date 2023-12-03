#!/usr/bin/env nextflow

// Using DSL-2
nextflow.enable.dsl=2

process train_torch_model {
    container "${params.container__mosdef_gomc}"
    publishDir "${params.output_folder}/torch", mode: 'copy', overwrite: true

    debug false
    input:
    file density_temperature_data
    output:
    tuple path("trained_model.pth"), path("scalers.joblib"), emit: model_scalers_tuple
    tuple path("predictions.png"), path("loss_epochs.png"), path("log.txt"), emit: model_performance
    script:
    """
    #!/usr/bin/env python
    # for python3
    import sys
    with open("log.txt", 'w') as sys.stdout:
        import pandas as pd
        import numpy as np
        import matplotlib.pyplot as plt
        from sklearn.model_selection import train_test_split
        from sklearn.preprocessing import StandardScaler
        from sklearn.metrics import mean_squared_error
        import torch
        import torch.nn as nn
        import torch.optim as optim
        from torch.utils.data import DataLoader, TensorDataset
        from torch.optim.lr_scheduler import ReduceLROnPlateau
        import joblib  # Import joblib for saving and loading scalers

        def save_scalers(scaler_X, scaler_y, scaler_file="scalers.joblib"):
            joblib.dump((scaler_X, scaler_y), scaler_file)
            print(f"Scalers saved to {scaler_file}")

        # Set seed for reproducibility
        seed = 42
        np.random.seed(seed)
        torch.manual_seed(seed)

        # Step 1: Read the CSV file into a pandas DataFrame
        file_path = "${density_temperature_data}"
        df = pd.read_csv(file_path)

        # Step 2: Extract the data for the neural network
        X = df["Density(kg/m3)"].values.reshape(-1, 1).astype(np.float32)
        y = df["Temperature(K)"].values.reshape(-1, 1).astype(np.float32)

        # Step 3: Split the data into training and testing sets
        X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=seed)

        # Step 4: Standardize the data
        scaler_X = StandardScaler()
        scaler_y = StandardScaler()

        X_train_scaled = scaler_X.fit_transform(X_train)
        X_test_scaled = scaler_X.transform(X_test)

        y_train_scaled = scaler_y.fit_transform(y_train)
        y_test_scaled = scaler_y.transform(y_test)

        # Step 5: Convert data to PyTorch tensors
        X_train_tensor = torch.from_numpy(X_train_scaled)
        y_train_tensor = torch.from_numpy(y_train_scaled)
        X_test_tensor = torch.from_numpy(X_test_scaled)
        y_test_tensor = torch.from_numpy(y_test_scaled)

        # Step 6: Create DataLoader for training
        train_dataset = TensorDataset(X_train_tensor, y_train_tensor)
        train_loader = DataLoader(train_dataset, batch_size=64, shuffle=True)

        # Step 7: Define the neural network model with more complexity and dropout
        class RegressionModel(nn.Module):
            def __init__(self, input_size):
                super(RegressionModel, self).__init__()
                self.fc1 = nn.Linear(input_size, 256)
                self.relu1 = nn.ReLU()
                self.dropout1 = nn.Dropout(0.2)  # Introducing dropout
                self.fc2 = nn.Linear(256, 128)
                self.relu2 = nn.ReLU()
                self.dropout2 = nn.Dropout(0.2)
                self.fc3 = nn.Linear(128, 64)
                self.relu3 = nn.ReLU()
                self.dropout3 = nn.Dropout(0.2)
                self.fc4 = nn.Linear(64, 1)

            def forward(self, x):
                x = self.fc1(x)
                x = self.relu1(x)
                x = self.dropout1(x)
                x = self.fc2(x)
                x = self.relu2(x)
                x = self.dropout2(x)
                x = self.fc3(x)
                x = self.relu3(x)
                x = self.dropout3(x)
                x = self.fc4(x)
                return x

        # Step 8: Instantiate the model, loss function, and optimizer with adjusted learning rate
        model = RegressionModel(input_size=1)
        criterion = nn.MSELoss()
        optimizer = optim.Adam(model.parameters(), lr=0.0005)  # Adjusted learning rate

        # Step 9: Learning rate scheduler
        scheduler = ReduceLROnPlateau(optimizer, mode='min', factor=0.1, patience=5, verbose=True)

        # Step 10: Train the model for more epochs with loss tracking
        num_epochs = 200
        losses = []

        for epoch in range(num_epochs):
            epoch_loss = 0.0
            for inputs, targets in train_loader:
                optimizer.zero_grad()
                outputs = model(inputs)
                loss = criterion(outputs, targets)
                loss.backward()
                optimizer.step()
                epoch_loss += loss.item()

            epoch_loss /= len(train_loader)
            losses.append(epoch_loss)

            if (epoch + 1) % 10 == 0:
                print(f"Epoch [{epoch + 1}/{num_epochs}], Loss: {epoch_loss:.4f}")

            # Update the learning rate based on the validation loss
            with torch.no_grad():
                model.eval()
                val_loss = criterion(model(X_test_tensor), y_test_tensor)
                scheduler.step(val_loss)

        # Step 11: Predict the full domain and plot it as a line
        X_domain = np.linspace(X.min(), X.max(), 100).reshape(-1, 1).astype(np.float32)
        X_domain_scaled = scaler_X.transform(X_domain)
        X_domain_tensor = torch.from_numpy(X_domain_scaled)
        with torch.no_grad():
            model.eval()
            y_domain_scaled = model(X_domain_tensor)
            y_domain = scaler_y.inverse_transform(y_domain_scaled.numpy())

        # Step 12: Evaluate the performance of the model on the test set
        with torch.no_grad():
            model.eval()
            y_pred_scaled = model(X_test_tensor)
            y_pred = scaler_y.inverse_transform(y_pred_scaled.numpy())
            mse = mean_squared_error(y_test, y_pred)
            print(f"Mean Squared Error on Test Set: {mse}")

        # Step 13: Plot the data and predictions
        plt.figure(figsize=(10, 6))
        plt.scatter(X, y, label="NIST miniREFPROP data", marker='o')
        plt.plot(X_domain, y_domain, color='red', label="Model Prediction", linewidth=2.5)  # Adjust linewidth
        plt.xlabel("Density (kg/m³)", fontsize=20)  # Doubled fontsize
        plt.ylabel("Temperature (K)", fontsize=20)  # Doubled fontsize
        plt.title("Water Density at 1 atm - Model Prediction and Data", fontsize=20)  # Doubled fontsize
        plt.legend(fontsize=16)  # Doubled fontsize
        plt.xticks(fontsize=16)  # Increased tick label size
        plt.yticks(fontsize=16)  # Increased tick label size
        plt.grid(True)
        plt.show()
        plt.savefig("predictions.png")

        # Step 14: Plot the training loss over epochs
        plt.figure(figsize=(10, 6))
        plt.plot(losses, label="Training Loss")
        plt.xlabel("Epoch", fontsize=20)
        plt.ylabel("Loss", fontsize=20)
        plt.title("Training Loss Over Epochs", fontsize=20)
        plt.legend(fontsize=16)
        plt.xticks(fontsize=16)  # Increased tick label size
        plt.yticks(fontsize=16)  # Increased tick label size
        plt.grid(True)
        plt.show()
        plt.savefig("loss_epochs.png")

        # Step 14: Save the trained model
        model_path = "trained_model.pth"
        torch.save(model.state_dict(), model_path)
        print(f"Model saved at {model_path}")

        # Instantiate and fit the scalers on training data if not already loaded
        scaler_file = "scalers.joblib"
        save_scalers(scaler_X, scaler_y, scaler_file)

    """
}


process predict_torch_model {
    container "${params.container__mosdef_gomc}"
    publishDir "${params.output_folder}/density_statepoints/${density}", mode: 'copy', overwrite: true

    debug false
    input:
    tuple path(torch_model), path(torch_scalers), val(density)
    output:
    tuple val(density), path("statepoint.json"), emit: statepoint

    script:
    """
    #!/usr/bin/env python
    import torch
    import numpy as np
    from sklearn.preprocessing import StandardScaler
    import joblib  # Import joblib for saving and loading scalers
    def load_scalers(scaler_file="scalers.joblib"):
        scaler_X, scaler_y = joblib.load(scaler_file)
        print(f"Scalers loaded from {scaler_file}")
        return scaler_X, scaler_y

    class RegressionModel(torch.nn.Module):
        def __init__(self, input_size):
            super(RegressionModel, self).__init__()
            self.fc1 = torch.nn.Linear(input_size, 256)
            self.relu1 = torch.nn.ReLU()
            self.dropout1 = torch.nn.Dropout(0.2)
            self.fc2 = torch.nn.Linear(256, 128)
            self.relu2 = torch.nn.ReLU()
            self.dropout2 = torch.nn.Dropout(0.2)
            self.fc3 = torch.nn.Linear(128, 64)
            self.relu3 = torch.nn.ReLU()
            self.dropout3 = torch.nn.Dropout(0.2)
            self.fc4 = torch.nn.Linear(64, 1)

        def forward(self, x):
            x = self.fc1(x)
            x = self.relu1(x)
            x = self.dropout1(x)
            x = self.fc2(x)
            x = self.relu2(x)
            x = self.dropout2(x)
            x = self.fc3(x)
            x = self.relu3(x)
            x = self.dropout3(x)
            x = self.fc4(x)
            return x

    def load_model(model_path, input_size=1):
        model = RegressionModel(input_size)
        model.load_state_dict(torch.load(model_path))
        model.eval()
        return model

    def predict_temperature(model, density, scaler_X, scaler_y):
        density_scaled = scaler_X.transform(np.array([[density]]).astype(np.float32))
        density_tensor = torch.from_numpy(density_scaled)
        
        with torch.no_grad():
            model.eval()
            temperature_scaled = model(density_tensor)
            temperature = scaler_y.inverse_transform(temperature_scaled.numpy())[0, 0]
        
        return temperature

    # Example usage:
    loaded_model = load_model("${torch_model}", input_size=1)
    scaler_X, scaler_y = load_scalers("${torch_scalers}")
    density_to_predict = ${density}  # Replace with the desired density
    predicted_temperature = predict_temperature(loaded_model, density_to_predict, scaler_X, scaler_y)
    print(f"Predicted Temperature for Density {density_to_predict}: {predicted_temperature:.2f} K")

    from typing import List
    from pydantic import BaseModel

    class Point(BaseModel):
        density: float
        temperature: float

    # Create a Pydantic object
    point_obj = Point(density=${density}, temperature=predicted_temperature)

    # Serialize the Pydantic object to JSON
    with open("statepoint.json", 'w') as file:
        file.write(point_obj.json())


    """
}


workflow train_model {
    take:
    csv_channel
    main:
    train_torch_model(csv_channel)
    emit:
    model_scalers_tuple = train_torch_model.out.model_scalers_tuple
    
}

workflow predict_model {
    take:
    model_density_tuple
    main:
    predict_torch_model(model_density_tuple)
    emit:
    statepoints = predict_torch_model.out.statepoint
    
}