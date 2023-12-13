import pandas as pd
from datetime import datetime, timedelta
import os

# This code takes 100 simulated call schedules as an input and returns an excel file where, for each simulation,
# an optimized schedule has been made. Furthermore, it returns the key parameters such as
# average waiting time and idle time for both patient/machine types in both the old and new system.
# Written by Joshua Smilde, Maciej Korek and Ties Bos

# Parameters for the scan durations
duration_type1 = 0.5 # MRI scan duration for Type 1 patients in hours 0.5
duration_type2 = 1  # MRI scan duration for Type 2 patients in hours 1.0

# Define operational hours
operation_start = 8  # 08:00 AM
operation_end = 17   # 05:00 PM

# Reading the CSV file into a DataFrame
df = pd.read_csv(r'Data/ScanRecords1.csv')

# Converting 'Date' and 'Time' to a single datetime object named 'CallDateTime'
df['CallDateTime'] = pd.to_datetime(df['Date']) + pd.to_timedelta(df['Time'], unit='h')

# Sorting the DataFrame by 'CallDateTime' to process in the order calls were received
df = df.sort_values('CallDateTime')

# Initialize the availability datetime for both MRI machines
mri1_availability = datetime.min
mri2_availability = datetime.min

# Lists to store scheduled scan times
scheduled_type1 = []
scheduled_type2 = []

def round_up_to_nearest_quarter_hour(dt):
    minutes = dt.minute
    if minutes % 15 != 0:
        dt = dt + timedelta(minutes=(15 - minutes % 15))
    return dt.replace(second=0, microsecond=0)

# Modified get_next_available_slot function
def get_next_available_slot(call_time, machine_availability, scan_duration):
    # Calculate the earliest possible scan time, which is the next day
    earliest_scan_time = (call_time + timedelta(days=1)).replace(hour=operation_start, minute=0, second=0, microsecond=0)
    # If the machine is available before the operation hours, set it to operation start time
    if machine_availability < earliest_scan_time:
        machine_availability = earliest_scan_time
    # Round up machine availability to the nearest quarter-hour
    machine_availability = round_up_to_nearest_quarter_hour(machine_availability)
    # If the scan fits in the operational hours of the day, schedule it
    if machine_availability.time() < datetime.min.replace(hour=operation_end).time():
        return max(earliest_scan_time, machine_availability)
    # Otherwise, schedule it for the next day's operation start time
    return round_up_to_nearest_quarter_hour(machine_availability.replace(hour=operation_start, minute=0, second=0, microsecond=0) + timedelta(days=1))

def schedule_scans(df, system_type, duration_type1, duration_type2):
    mri1_availability = datetime.min
    mri2_availability = datetime.min
    scheduled_type1 = []
    scheduled_type2 = []
    scheduled_mri1 = []
    scheduled_mri2 = []

    for index, row in df.iterrows():
        call_time = row['CallDateTime']
        if system_type == 'old':
            if row['PatientType'] == 'Type 1':
                # Schedule for MRI1 in Old System
                scheduled_time = get_next_available_slot(call_time, mri1_availability, duration_type1)
                scheduled_type1.append((row['CallDateTime'], row['PatientType'], scheduled_time))
                mri1_availability = scheduled_time + timedelta(hours=duration_type1)
            else:
                # Schedule for MRI2 in Old System
                scheduled_time = get_next_available_slot(call_time, mri2_availability, duration_type2)
                scheduled_type2.append((row['CallDateTime'], row['PatientType'], scheduled_time))
                mri2_availability = scheduled_time + timedelta(hours=duration_type2)
        else:
            # New system scheduling logic
            scan_duration = duration_type1 if row['PatientType'] == 'Type 1' else duration_type2
            if mri1_availability <= mri2_availability:
                scheduled_time = get_next_available_slot(call_time, mri1_availability, scan_duration)
                scheduled_mri1.append((row['CallDateTime'], row['PatientType'], scheduled_time))
                mri1_availability = scheduled_time + timedelta(hours=scan_duration)
            else:
                scheduled_time = get_next_available_slot(call_time, mri2_availability, scan_duration)
                scheduled_mri2.append((row['CallDateTime'], row['PatientType'], scheduled_time))
                mri2_availability = scheduled_time + timedelta(hours=scan_duration)

    if system_type == 'old':
        return pd.DataFrame(scheduled_type1, columns=['CallDateTime', 'PatientType', 'ScheduledTime']), pd.DataFrame(scheduled_type2, columns=['CallDateTime', 'PatientType', 'ScheduledTime'])
    else:
        return pd.DataFrame(scheduled_mri1, columns=['CallDateTime', 'PatientType', 'ScheduledTime']), pd.DataFrame(scheduled_mri2, columns=['CallDateTime', 'PatientType', 'ScheduledTime'])

def calculate_waiting_time(call_time, scheduled_time):
    waiting_time = 0

    while call_time < scheduled_time:
        # Define next_hour at the beginning of the loop
        next_hour = call_time.replace(minute=0, second=0, microsecond=0) + timedelta(hours=1)

        # If within operational hours, calculate partial hour if needed
        if operation_start <= call_time.hour < operation_end:
            if scheduled_time < next_hour:
                waiting_time += (scheduled_time - call_time).total_seconds() / 3600
                break
            else:
                waiting_time += (next_hour - call_time).total_seconds() / 3600

        # Move to the next hour or next operational day
        call_time = next_hour
        if call_time.hour >= operation_end:
            call_time = call_time.replace(hour=operation_start, minute=0, second=0, microsecond=0) + timedelta(days=1)

    return waiting_time

def calculate_waiting_times(scheduled_dfs):
    waiting_times = []

    for scheduled_df in scheduled_dfs:
        for _, row in scheduled_df.iterrows():
            call_time = row['CallDateTime']
            scheduled_time = row['ScheduledTime']
            waiting_time = calculate_waiting_time(call_time, scheduled_time)
            waiting_times.append(waiting_time)

    avg_waiting_time = sum(waiting_times) / len(waiting_times) if waiting_times else 0
    max_waiting_time = max(waiting_times) if waiting_times else 0

    return avg_waiting_time, max_waiting_time

def calculate_idle_time(scheduled_df):
    last_scan_end = datetime.min
    total_idle_time_seconds = 0

    for _, row in scheduled_df.iterrows():
        scheduled_start = row['ScheduledTime']
        if last_scan_end != datetime.min and operation_start <= last_scan_end.hour < operation_end:
            idle_time = scheduled_start - last_scan_end
            total_idle_time_seconds += max(idle_time.total_seconds(), 0)
        last_scan_end = scheduled_start + timedelta(hours=row['ScanDuration'])

    operational_days = (scheduled_df['ScheduledTime'].dt.date.max() - scheduled_df['ScheduledTime'].dt.date.min()).days + 1
    avg_idle_time_hours = (total_idle_time_seconds / 3600) / operational_days if operational_days > 0 else 0

    return avg_idle_time_hours


def calculate_average_appointments_per_day(scheduled_dfs):
    combined_df = pd.concat(scheduled_dfs)
    combined_df['Date'] = combined_df['ScheduledTime'].dt.date

    daily_counts = combined_df.groupby(['Date', 'PatientType']).size().reset_index(name='Counts')
    average_counts = daily_counts.groupby('PatientType')['Counts'].mean()

    return average_counts

def calculate_total_appointments(df):
    total_counts = df['PatientType'].value_counts()
    return total_counts


def analyze_scan_records(df):
    # Convert 'Date' and 'Time' to datetime object 'CallDateTime'
    df['CallDateTime'] = pd.to_datetime(df['Date']) + pd.to_timedelta(df['Time'], unit='h')
    df = df.sort_values('CallDateTime')

    # Old System Analysis
    scheduled_old_type1, scheduled_old_type2 = schedule_scans(df, 'old', duration_type1, duration_type2)
    scheduled_old_type1['ScanDuration'] = duration_type1
    scheduled_old_type2['ScanDuration'] = duration_type2

    # New System Analysis
    scheduled_new_mri1, scheduled_new_mri2 = schedule_scans(df, 'new', duration_type1, duration_type2)
    scheduled_new_mri1['ScanDuration'] = scheduled_new_mri1['PatientType'].apply(
        lambda x: duration_type1 if x == 'Type 1' else duration_type2)
    scheduled_new_mri2['ScanDuration'] = scheduled_new_mri2['PatientType'].apply(
        lambda x: duration_type1 if x == 'Type 1' else duration_type2)

    # Calculate waiting times for each patient type in the Old System
    waiting_times_old_type1 = [calculate_waiting_time(row['CallDateTime'], row['ScheduledTime']) for index, row in
                               scheduled_old_type1.iterrows()]
    waiting_times_old_type2 = [calculate_waiting_time(row['CallDateTime'], row['ScheduledTime']) for index, row in
                               scheduled_old_type2.iterrows()]

    # Calculate waiting times for each patient type in the New System
    waiting_times_new_type1 = [calculate_waiting_time(row['CallDateTime'], row['ScheduledTime']) for index, row in
                               scheduled_new_mri1.iterrows()]
    waiting_times_new_type2 = [calculate_waiting_time(row['CallDateTime'], row['ScheduledTime']) for index, row in
                               scheduled_new_mri2.iterrows()]

    # Function to calculate average and maximum waiting times
    def calculate_avg_max(waiting_times):
        avg_waiting = sum(waiting_times) / len(waiting_times) if waiting_times else 0
        max_waiting = max(waiting_times) if waiting_times else 0
        return avg_waiting, max_waiting

    # Calculating average and maximum waiting times
    avg_waiting_old_type1, max_waiting_old_type1 = calculate_avg_max(waiting_times_old_type1)
    avg_waiting_old_type2, max_waiting_old_type2 = calculate_avg_max(waiting_times_old_type2)
    avg_waiting_new_type1, max_waiting_new_type1 = calculate_avg_max(waiting_times_new_type1)
    avg_waiting_new_type2, max_waiting_new_type2 = calculate_avg_max(waiting_times_new_type2)

    # Calculate Idle Times separately for each machine
    idle_time_old_mri1 = calculate_idle_time(scheduled_old_type1)
    idle_time_old_mri2 = calculate_idle_time(scheduled_old_type2)
    idle_time_new_mri1 = calculate_idle_time(scheduled_new_mri1)
    idle_time_new_mri2 = calculate_idle_time(scheduled_new_mri2)

    # Calculate Average Appointments per Day
    avg_appointments = calculate_average_appointments_per_day(
        [scheduled_old_type1, scheduled_old_type2, scheduled_new_mri1, scheduled_new_mri2])

    # Calculate Total Appointments
    total_appointments = calculate_total_appointments(df)

    # Print all the results
    print("\nOld System Waiting Times:")
    print(f"Type 1 - Average: {avg_waiting_old_type1}, Max: {max_waiting_old_type1}")
    print(f"Type 2 - Average: {avg_waiting_old_type2}, Max: {max_waiting_old_type2}")

    print("\nNew System Waiting Times:")
    print(f"Type 1 - Average: {avg_waiting_new_type1}, Max: {max_waiting_new_type1}")
    print(f"Type 2 - Average: {avg_waiting_new_type2}, Max: {max_waiting_new_type2}")

    print("\nIdle Time - Old System MRI 1:", idle_time_old_mri1)
    print("Idle Time - Old System MRI 2:", idle_time_old_mri2)
    print("Idle Time - New System MRI 1:", idle_time_new_mri1)
    print("Idle Time - New System MRI 2:", idle_time_new_mri2)

    print("\nAverage Appointments per Day:")
    print(avg_appointments)

    print("\nTotal Appointments for Each Patient Type:")
    print(total_appointments)

    return {
        'avg_waiting_old_type1': avg_waiting_old_type1,
        'max_waiting_old_type1': max_waiting_old_type1,
        'avg_waiting_old_type2': avg_waiting_old_type2,
        'max_waiting_old_type2': max_waiting_old_type2,
        'avg_waiting_new_type1': avg_waiting_new_type1,
        'max_waiting_new_type1': max_waiting_new_type1,
        'avg_waiting_new_type2': avg_waiting_new_type2,
        'max_waiting_new_type2': max_waiting_new_type2,
        'idle_time_old_mri1': idle_time_old_mri1,
        'idle_time_old_mri2': idle_time_old_mri2,
        'idle_time_new_mri1': idle_time_new_mri1,
        'idle_time_new_mri2': idle_time_new_mri2,
        'avg_appointments': avg_appointments,
        'total_appointments': total_appointments
    }

results = analyze_scan_records(df)

results_list = []

for i in range(1, 101):  # Assuming you have files from ScanRecords1 to ScanRecords100
    file_path = f'Data/ScanRecords{i}.csv'
    if os.path.exists(file_path):
        df = pd.read_csv(file_path)
        df['CallDateTime'] = pd.to_datetime(df['Date']) + pd.to_timedelta(df['Time'], unit='h')
        df = df.sort_values('CallDateTime')

        results = analyze_scan_records(df)
        results['file'] = f'ScanRecords{i}'
        results_list.append(results)

# Convert the list of results to a DataFrame
results_df = pd.DataFrame(results_list)

# Export to Excel
results_df.to_excel('ScanRecords_Analysis.xlsx', index=False)