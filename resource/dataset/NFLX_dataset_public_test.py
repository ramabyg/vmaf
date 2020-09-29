dataset_name = 'NFLX_public_additnal_feature'
yuv_fmt = 'yuv420p'
width = 1920
height = 1080
eotf = 1
clrFmt = 2

ref_dir = '/home/rama/rama/vmaf/data/NFLX_dataset/ref'
dis_dir = '/home/rama/rama/vmaf/data/NFLX_dataset/dis'

ref_videos = [
 {'content_id': 0,
  'content_name': 'BigBuckBunny',
  'path': ref_dir + '/BigBuckBunny_25fps.yuv'},
 {'content_id': 1,
  'content_name': 'BirdsInCage',
  'path': ref_dir + '/BirdsInCage_30fps.yuv'},
 {'content_id': 2,
  'content_name': 'CrowdRun',
  'path': ref_dir + '/CrowdRun_25fps.yuv'},
 {'content_id': 3,
  'content_name': 'ElFuente1',
  'path': ref_dir + '/ElFuente1_30fps.yuv'},
 {'content_id': 4,
  'content_name': 'ElFuente2',
  'path': ref_dir + '/ElFuente2_30fps.yuv'},
 {'content_id': 5,
  'content_name': 'FoxBird',
  'path': ref_dir + '/FoxBird_25fps.yuv'},
 {'content_id': 6,
  'content_name': 'OldTownCross',
  'path': ref_dir + '/OldTownCross_25fps.yuv'},
 {'content_id': 7,
  'content_name': 'Seeking',
  'path': ref_dir + '/Seeking_25fps.yuv'},
 {'content_id': 8,
  'content_name': 'Tennis',
  'path': ref_dir + '/Tennis_24fps.yuv'}
]

dis_videos = [
 {'asset_id': 14,
  'content_id': 0,
  'dmos': 97.5,
  'path': dis_dir + '/BigBuckBunny_75_720_3050.yuv',
 }

]
