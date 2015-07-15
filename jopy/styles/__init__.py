
if __name__ == "__main__":
    from jopy.styles.mplib import IpuStyle
    line_fig,map_fig = IpuStyle()._show_info()
    line_fig.savefig("IpuStyle_lines.pdf")
    map_fig.savefig("IpuStyle_maps.pdf")