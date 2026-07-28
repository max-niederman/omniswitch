-- Trivial headless-run check: file IO works and FEMM exits cleanly.
-- FEMM embeds Lua 4.0 with all stdlib functions in the global namespace.
handle = openfile("Z:/home/max/Projects/hw/omniswitch/results/hello.txt", "w")
write(handle, "ok\n")
closefile(handle)
quit()
