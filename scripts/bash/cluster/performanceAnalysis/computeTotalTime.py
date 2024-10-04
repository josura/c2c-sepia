import pandas as pd

# Read the data, the data has the following columns:
# job_id	user	partition	date	duration	nodes
# 7117930	glocicer	boost_usr+	2024-08-20T17:32:44	00:00:03	4
# 7168517	glocicer	boost_usr+	2024-08-20T21:32:38	00:00:02	4
# 7168518	glocicer	boost_usr+	2024-08-20T21:32:38	00:00:02	4
# we are interested in nodes and duration, that should be converted to hours

data = pd.read_csv('/home/josura/Projects/ccc/c2c-sepia/scripts/bash/cluster/performanceAnalysis/jobsSubmitted.tsv', sep='\t')
data['duration'] = pd.to_timedelta(data['duration'])
data['duration'] = data['duration'].dt.total_seconds() / 3600
data['total_time'] = data['duration'] * data['nodes'] * 32
print(data['total_time'].sum())
