from qiskit import register, available_backends , get_backend 

# Establish connection with IBMQuantum Experience 
try :  
	import sys 
	sys.path.append('../')
	import Qconfig 
	qx_config = {                              # configuration details 
		'APItoken' : Qconfig.APItoken ,
 		'url' : Qconfig.config['url']}
except Exception as e :
	print(e)
	print("Check your API token") 

register(qx_config['APItoken'],qx_config['url'])

def lowest_pending_jobs(): # find the best backend available 
	list_backends = available_backends({'local':False,'simulator':False})
	device_status = [get_backend(backend).status for backend in list_backends]
	best = min([x for x in device_status if x['available'] is True],key = lambda x: x['pending_jobs'])
	return best['name']

def get_qc():
	backend = lowest_pending_jobs() 
	print("The best backend is",backend)
	return backend

 





