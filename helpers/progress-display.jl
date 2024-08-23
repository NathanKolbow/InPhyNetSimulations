function print_monophyleticRobustnessProgress(iterspassed::Int64, totaliters::Int64, start_time::Float64)
    percent_complete = round(100 * iterspassed / totaliters, digits = 2)
    time_passed = time() - start_time
    est_total_runtime = time_passed / (percent_complete / 100)
    units = "seconds"

    # if more than 90 min passed, print in hours
    # elseif more than 180 sec passed, print in minutes
    if est_total_runtime > 60 * 90
        time_passed /= 60 * 60
        est_total_runtime /= 60 * 60
        units = "hours"
    elseif est_total_runtime > 180
        time_passed /= 60
        est_total_runtime /= 60
        units = "minutes"
    end
    time_passed = round(time_passed, digits = 2)
    est_total_runtime = round(est_total_runtime, digits = 2)

    print("\r\t$(percent_complete)% ($(iterspassed)/$(totaliters)) complete "*
          "($(time_passed)/≈$(est_total_runtime) $(units))              ")
end


function display_progress(iterspassed::Int64, totaliters::Int64, start_time::Float64; npegs::Int64=70)
    if iterspassed == 0
        prog_header = "⊻ PROGRESS ⊻"
        print("\n\t ")
        for j = 1:((npegs ÷ 2) - (length(prog_header) ÷ 2)) print(" ") end
        print("$(prog_header)\n")

        print("0%\t| ")
        for j = 1:npegs print("-") end
        println(" |\t100%")
    else
        print("\r\t| ")

        if iterspassed == totaliters
            for j = 1:npegs print("*") end
        else
            n_draw = (npegs * iterspassed) ÷ totaliters
            for j = 1:n_draw print("*") end
            for j = 1:(npegs - n_draw) print(" ") end
        end
        
        print(" |\t$(iterspassed)/$(totaliters)  ")

        if iterspassed == totaliters
            # Print total time taken
            time_seconds = time() - start_time
            time_minutes = time_seconds / 60
            time_hours = time_minutes / 60
            time_days = time_hours / 24

            time_footer = ""
            if time_minutes < 1
                time_val = round(time_seconds, digits = 2)
                time_footer = string(time_val)
                time_footer = "$(time_footer) seconds"
            elseif time_hours < 1
                time_val = round(time_minutes, digits = 2)
                time_footer = string(time_val)
                time_footer = "$(time_footer) minutes"
            elseif time_days < 1
                time_val = round(time_hours, digits = 2)
                time_footer = string(time_val)
                time_footer = "$(time_footer) hours"
            else
                time_val = round(time_days, digits = 2)
                time_footer = string(time_val)
                time_footer = "$(time_footer) days"
            end
            
            print("\n\t ")
            for j = 1:((npegs ÷ 2) - (length(time_footer) ÷ 2)) print(" ") end
            print("$(time_footer)\n")

            
            # Print average time taken
            time_seconds /= totaliters
            time_minutes /= totaliters
            time_hours /= totaliters
            time_days /= totaliters

            time_footer = ""
            if time_minutes < 1
                time_val = round(time_seconds, digits = 2)
                time_footer = string(time_val)
                time_footer = "μ = $(time_footer) seconds"
            elseif time_hours < 1
                time_val = round(time_minutes, digits = 2)
                time_footer = string(time_val)
                time_footer = "μ = $(time_footer) minutes"
            elseif time_days < 1
                time_val = round(time_hours, digits = 2)
                time_footer = string(time_val)
                time_footer = "μ = $(time_footer) hours"
            else
                time_val = round(time_days, digits = 2)
                time_footer = string(time_val)
                time_footer = "μ = $(time_footer) days"
            end
            
            print("\t ")
            for j = 1:((npegs ÷ 2) - (length(time_footer) ÷ 2)) print(" ") end
            print("$(time_footer)\n\n")
        end
    end
end