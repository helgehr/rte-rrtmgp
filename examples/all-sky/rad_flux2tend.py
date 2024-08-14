import xarray as xr
import sys

assert len(sys.argv) == 4, "Need <input_sw> <input_lw> <output>"

cpd   = 1004.64 # J K-1 kg-1
g = 9.80616 # m s-2

ds_sw = xr.open_dataset(sys.argv[1])
ds_lw = xr.open_dataset(sys.argv[2])

out_ds = xr.Dataset()

net_rad_sw = ds_sw.sw_flux_up - ds_sw.sw_flux_dn
net_rad_lw = ds_lw.lw_flux_up - ds_lw.lw_flux_dn
net_rad = net_rad_sw + net_rad_lw
heating_rate = net_rad.diff(dim='lev') / ds_sw.p_lev.diff(dim='lev') * g / cpd

heating_rate = heating_rate.rename('ptend_t')

# climsim out
out_ds['ptend_t'] = heating_rate
out_ds['NETSW'] = -net_rad_sw.isel(lev=-1) # apparently defined as negative netflux (sfc absorption)
out_ds['FLWDS'] = ds_lw.lw_flux_dn.isel(lev=-1)
out_ds['SOLS'] = ds_sw.sw_flux_dir.isel(lev=-1)
# out_ds['SOLL']
# out_ds['SOLSD']
# out_ds['SOLLD']

# climsim in
out_ds['LWUP'] = ds_lw.lw_flux_up.isel(lev=-1)

# fix_units
out_ds['ptend_t'].attrs['Units'] = 'K s-1'
out_ds['NETSW'].attrs['Units'] = 'W m-2'
out_ds['FLWDS'].attrs['Units'] = 'W m-2'
out_ds['SOLS'].attrs['Units'] = 'W m-2'
out_ds['LWUP'].attrs['Units'] = 'W m-2'


out_ds.to_netcdf(sys.argv[3])
print(f'Wrote radiation output to {sys.argv[3]}')

# heating_rate.to_netcdf(sys.argv[3])
