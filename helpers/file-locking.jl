using Pidfile

function lock_file(out_file::String)
    lk = mkpidlock(out_file*".lk")
end