nbr_chunks = 8
samples_per_chunk = 180
number_channels = 2
tau_prc = 0.1
frame_rate = 4.
tau = 1/frame_rate - tau_prc
total_samples = nbr_chunks*samples_per_chunk*number_channels
Fs = total_samples/tau
RPM = 60*nbr_chunks/tau
print(f'Want {nbr_chunks} chunks, each having {samples_per_chunk} samples in {1000*tau} ms.')
print(f'Need to collect {total_samples/number_channels} samples per channel for a total of {total_samples} samples.')
print(f'Collect at {Fs/1000} ksps. Rotate at {RPM} RPM')