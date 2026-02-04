
plt.figure()
for clamp, style in thf_conditions:
    t_s, ch2 = simulate_ch2thf_with_thf_clamp_zigzag(clamp_thf=clamp)
    # keep all curves black but distinguish by pattern/markers
    kwargs = dict(color="black", linewidth=1.2, zorder=2)
    if style["marker"]:
        kwargs.update(marker=style["marker"], markevery=60, markersize=3)
    plt.plot(t_s, ch2, linestyle=style["linestyle"], label=style["label"], **kwargs)

# reference line at zero to highlight negative CH2THF if/when it occurs
plt.axhline(0, linewidth=0.8, color="black", zorder=1)

plt.xlabel("Time (s)")
plt.ylabel("CH2THF (µM)")
plt.title("CH2THF vs Time for Different THF Clamps (THF forced each step)")
plt.legend()
plt.tight_layout()
plt.show()