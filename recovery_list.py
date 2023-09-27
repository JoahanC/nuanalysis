with open("binning_check_3.txt", 'r') as file:
    with open("recovery_list_3.txt", 'a') as recovery:
        for line in file.readlines():
            data = line.split()
            chan_1 = float(data[3])
            chan_2 = float(data[4])
            chan_3 = float(data[5])
            chan_4 = float(data[6])
            chans = [chan_1, chan_2, chan_3, chan_4]
            flag = False
            for chan in chans:
                if flag:
                    continue
                if chan < 1.0:
                    recovery.write(line)
                    flag = True

