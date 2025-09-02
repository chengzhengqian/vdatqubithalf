for norb in 2:5
    for r in [0,0.05,0.1,0.25]
        print("norb $(norb) r $(r)\n")
        tag="two_peak"
        saveData(solve_r_norb(r,norb,cal_Δ_two_peak,cal_A_two_peak,cal_K_two_peak,data_dir;ainit=1.9,binit=0.1),filename_results(r,norb,tag,data_dir))
    end
end
