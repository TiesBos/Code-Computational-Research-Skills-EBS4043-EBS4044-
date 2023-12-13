import numpy as np
import pandas as pd
from datetime import datetime, timedelta

# This code creates 100 schedules for calls coming in, overa  period of 23 days.
# It contains calls from both patient types, and takes the relevant distributions,
# for the amount of people calling per day and the average time between the calls, into account.
# Written by Joshua Smilde, Maciej Korek and Ties Bos. 

# Parameters for the simulation
mean_duration_type1 = 0.4327
std_dev_duration_type1 = 0.0978

lambda_type1 = 379 / 23
lambda_type2 = 239 / 23

# Gamma distribution parameters for Type 2
alpha_type2 = 12.77268
beta_type2 = 19.0825

# Function to generate weekdays in August
def generate_weekdays_in_august(year):
    start_date = datetime(year, 8, 1)
    end_date = datetime(year, 8, 31)
    return [start_date + timedelta(days=x) for x in range((end_date - start_date).days + 1) if (start_date + timedelta(days=x)).weekday() < 5]

# Generating weekdays for August 2023
weekdays_in_august = generate_weekdays_in_august(2023)

# Modify the simulate_interleaved_schedule function to match the format of the uploaded file
def simulate_interleaved_schedule_modified():
    records = []

    for date in weekdays_in_august:
        # Simulating number of calls for each type
        num_calls_type1 = np.random.poisson(lambda_type1)
        num_calls_type2 = np.random.poisson(lambda_type2)

        # Initial call time starts at 8:00
        call_time_type1 = 8.0
        call_time_type2 = 8.0

        # Generating call records for Type 1
        for _ in range(num_calls_type1):
            duration = max(0, np.random.normal(mean_duration_type1, std_dev_duration_type1))
            records.append([date.strftime('%Y-%m-%d'), call_time_type1, duration, 'Type 1'])
            # Update call time
            call_time_type1 += np.random.exponential(1 / 1.833705)

        # Generating call records for Type 2
        for _ in range(num_calls_type2):
            duration = max(0, np.random.gamma(alpha_type2, 1/beta_type2))
            records.append([date.strftime('%Y-%m-%d'), call_time_type2, duration, 'Type 2'])
            # Update call time
            call_time_type2 += max(0, np.random.normal(0.8666, 0.31078))

    # Combining and sorting Type 1 and Type 2 calls
    combined_schedule = pd.DataFrame(records, columns=['Date', 'Time', 'Duration', 'PatientType'])
    combined_schedule.sort_values(by=['Date', 'Time'], inplace=True)

    return combined_schedule

def create_and_save_dataset(file_number):
    # Generate the schedule
    interleaved_schedule = simulate_interleaved_schedule_modified()

    # Calculate and print the average number of calls per day for each patient type
    average_calls_per_day_type1 = interleaved_schedule[interleaved_schedule['PatientType'] == 'Type 1'].groupby('Date').size().mean()
    average_calls_per_day_type2 = interleaved_schedule[interleaved_schedule['PatientType'] == 'Type 2'].groupby('Date').size().mean()

    print(f"Average calls per day for Type 1 in ScanRecords{file_number}: {average_calls_per_day_type1}")
    print(f"Average calls per day for Type 2 in ScanRecords{file_number}: {average_calls_per_day_type2}")

    # Path for the new file
    output_file_path = rf'Data/ScanRecords{file_number}.csv'

    # Write to a CSV file
    interleaved_schedule.to_csv(output_file_path, index=False)

# Main loop to create and save 100 datasets
for i in range(1, 101):
    create_and_save_dataset(i)
