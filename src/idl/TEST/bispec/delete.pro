pro delete, filename
case !version.os of
    'vms': cmd = 'delete'
    'windows': cmd = 'del'
    'MacOS': goto, notsup
    else: cmd = '\rm -f'
    endcase
cmd = cmd + ' ' + filename  
spawn, cmd  
goto, DONE
NOTSUP: print, "This operation is not supported on the Macintosh'
DONE:
end


