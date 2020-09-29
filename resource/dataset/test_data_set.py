dataset_name = 'Test_Data_set'
yuv_fmt = 'yuv420p'
width = 576
height = 324
#eotf = 1
clrFmt = 2

ref_dir = '/home/rama/rama/vmaf/data/yuv'
dis_dir = '/home/rama/rama/vmaf/data/yuv'

ref_videos = [
 {'content_id': 0,
  'content_name': 'test1',
  'path': ref_dir + '/src01_hrc00_576x324.yuv'},
 {'content_id': 1,
  'content_name': 'test2',
  'path': ref_dir + '/src01_hrc00_576x324.yuv'},
 {'content_id': 2,
  'content_name': 'test3',
  'path': ref_dir + '/src01_hrc00_576x324.yuv'},
]

dis_videos = [
 {'asset_id': 1,
  'content_id': 0,
  'dmos': 97.5,
  'path': dis_dir + '/src01_hrc01_576x324.yuv',
 },
 
 {'asset_id': 2,
    'content_id': 1,
    'dmos': 97.5,
    'path': dis_dir + '/src01_hrc01_576x324.yuv',
 }
]
