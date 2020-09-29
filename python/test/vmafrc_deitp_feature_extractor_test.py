import unittest

from vmaf.core.asset import Asset

from vmaf.config import VmafConfig

from test.testutil import set_default_576_324_videos_for_testing, set_default_576_324_12bit_videos_for_testing, \
    set_default_576_324_16bit_videos_for_testing, set_default_576_324_10bit_videos_for_testing_b, \
    set_default_576_324_10bit_videos_for_testing, set_default_576_324_yuv422p_videos_for_testing, \
    set_default_576_324_yuv444p_videos_for_testing, set_default_576_324_yuv444p10le_videos_for_testing, \
    set_default_576_324_yuv444p12le_videos_for_testing

from vmaf.core.feature_extractor import FloatDeitpFeatureExtractor


class DeitpFeatureExtractorTest(unittest.TestCase):

    def tearDown(self):
        if hasattr(self, 'fextractor'):
            self.fextractor.remove_results()
        pass

    def test_float_deitp_yuv420_8bit(self):
        ref_path, dis_path, asset, asset_original = set_default_576_324_videos_for_testing()

        # Test with Full Range
        self.fextractor = FloatDeitpFeatureExtractor(
            [asset, asset_original],
            None, fifo_mode=False,
            result_store=None,
            optional_dict={'eotf': 0, 'gamma': 2.4, 'clrFmt': 0, 'minLum': 0.005,
                           'maxLum': 100.0, 'Yuv2RgbExt': 'false',
                           'YuvXferSpec': 2, 'Rng': 1, 'RgbDef': 2}
            #optional_dict={'eotf': 0, 'gamma': 2.4}
        )
        self.fextractor.run()
        results = self.fextractor.results

       # print(results[0])
       #  print(results[1])

        self.assertAlmostEqual(results[0]['DEITP_feature_deitp_max_score'], 119.089463, places=5)
        self.assertAlmostEqual(results[0]['DEITP_feature_deitp_mean_score'], 12.284851, places=5)
        self.assertAlmostEqual(results[0]['DEITP_feature_deitp_sd_score'], 7.693644, places=5)
        self.assertAlmostEqual(results[1]['DEITP_feature_deitp_max_score'], 0.000000, places=6)
        self.assertAlmostEqual(results[1]['DEITP_feature_deitp_mean_score'], 0.000000, places=6)
        self.assertAlmostEqual(results[1]['DEITP_feature_deitp_sd_score'], 0.000000, places=6)
        
        self.assertAlmostEqual(results[0]['DEITP_feature_deTP_max_score'], 68.137123, places=5)
        self.assertAlmostEqual(results[0]['DEITP_feature_deTP_mean_score'], 5.750864, places=5)
        self.assertAlmostEqual(results[0]['DEITP_feature_deTP_sd_score'], 3.522399, places=5)
        self.assertAlmostEqual(results[1]['DEITP_feature_deTP_max_score'], 0.000000, places=6)
        self.assertAlmostEqual(results[1]['DEITP_feature_deTP_mean_score'], 0.000000, places=6)
        self.assertAlmostEqual(results[1]['DEITP_feature_deTP_sd_score'], 0.000000, places=6)
    
        #with default parameters
        self.fextractor = FloatDeitpFeatureExtractor(
            [asset, asset_original],
            None, fifo_mode=False,
            result_store=None,
            optional_dict={'eotf': 0, 'gamma': 2.4, 'clrFmt': 0, 'minLum': 0.005, 
                            'maxLum': 100.0, 'Yuv2RgbExt': 'false',
                            'YuvXferSpec': 2, 'Rng': 0, 'RgbDef':2}
        )
        self.fextractor.run()
        results = self.fextractor.results

        #print(results[0])
      #  print(results[1])

        self.assertAlmostEqual(results[0]['DEITP_feature_deitp_max_score'], 145.079002, places=5)
        self.assertAlmostEqual(results[0]['DEITP_feature_deitp_mean_score'], 14.733384, places=5)
        self.assertAlmostEqual(results[0]['DEITP_feature_deitp_sd_score'], 9.864875, places=5)
        self.assertAlmostEqual(results[1]['DEITP_feature_deitp_max_score'], 0.000000, places=6)
        self.assertAlmostEqual(results[1]['DEITP_feature_deitp_mean_score'], 0.000000, places=6)
        self.assertAlmostEqual(results[1]['DEITP_feature_deitp_sd_score'], 0.000000, places=6)

        self.assertAlmostEqual(results[0]['DEITP_feature_deTP_max_score'], 79.211248, places=5)
        self.assertAlmostEqual(results[0]['DEITP_feature_deTP_mean_score'], 6.164932, places=5)
        self.assertAlmostEqual(results[0]['DEITP_feature_deTP_sd_score'], 4.017366, places=5)
        self.assertAlmostEqual(results[1]['DEITP_feature_deTP_max_score'], 0.000000, places=6)
        self.assertAlmostEqual(results[1]['DEITP_feature_deTP_mean_score'], 0.000000, places=6)
        self.assertAlmostEqual(results[1]['DEITP_feature_deTP_sd_score'], 0.000000, places=6)
        
    def test_float_deitp_yuv420_10bit(self):
        ref_path, dis_path, asset, asset_original = set_default_576_324_10bit_videos_for_testing_b()
        self.fextractor = FloatDeitpFeatureExtractor(
            [asset, asset_original],
            None, fifo_mode=False,
            result_store=None
        )
        self.fextractor.run()
        results = self.fextractor.results

       # print(results[0])
       # print(results[1])

        self.assertAlmostEqual(results[0]['DEITP_feature_deitp_max_score'], 117.734121, places=5)
        self.assertAlmostEqual(results[0]['DEITP_feature_deitp_mean_score'], 12.695032, places=5)
        self.assertAlmostEqual(results[0]['DEITP_feature_deitp_sd_score'], 8.085674, places=5)
        self.assertAlmostEqual(results[1]['DEITP_feature_deitp_max_score'], 0.000000, places=6)
        self.assertAlmostEqual(results[1]['DEITP_feature_deitp_mean_score'], 0.000000, places=6)
        self.assertAlmostEqual(results[1]['DEITP_feature_deitp_sd_score'], 0.000000, places=6)
        
        self.assertAlmostEqual(results[0]['DEITP_feature_deTP_max_score'], 59.954721, places=5)
        self.assertAlmostEqual(results[0]['DEITP_feature_deTP_mean_score'], 5.902296, places=5)
        self.assertAlmostEqual(results[0]['DEITP_feature_deTP_sd_score'], 3.678506, places=5)
        self.assertAlmostEqual(results[1]['DEITP_feature_deTP_max_score'], 0.000000, places=6)
        self.assertAlmostEqual(results[1]['DEITP_feature_deTP_mean_score'], 0.000000, places=6)
        self.assertAlmostEqual(results[1]['DEITP_feature_deTP_sd_score'], 0.000000, places=6)

    def test_float_deitp_yuv420_16bit(self):
        ref_path, dis_path, asset, asset_original = set_default_576_324_16bit_videos_for_testing()
        self.fextractor = FloatDeitpFeatureExtractor(
            [asset, asset_original],
            None, fifo_mode=False,
            result_store=None
        )
        self.fextractor.run()
        results = self.fextractor.results

        #print(results[0])
        # print(results[1])


        self.assertAlmostEqual(results[0]['DEITP_feature_deitp_max_score'], 117.734121, places=5)
        self.assertAlmostEqual(results[0]['DEITP_feature_deitp_mean_score'], 12.695051, places=5)
        self.assertAlmostEqual(results[0]['DEITP_feature_deitp_sd_score'], 8.085678, places=5)
        self.assertAlmostEqual(results[1]['DEITP_feature_deitp_max_score'], 0.000000, places=6)
        self.assertAlmostEqual(results[1]['DEITP_feature_deitp_mean_score'], 0.000000, places=6)
        self.assertAlmostEqual(results[1]['DEITP_feature_deitp_sd_score'], 0.000000, places=6)

        self.assertAlmostEqual(results[0]['DEITP_feature_deTP_max_score'], 59.955428, places=5)
        self.assertAlmostEqual(results[0]['DEITP_feature_deTP_mean_score'], 5.902295, places=5)
        self.assertAlmostEqual(results[0]['DEITP_feature_deTP_sd_score'], 3.678510, places=5)
        self.assertAlmostEqual(results[1]['DEITP_feature_deTP_max_score'], 0.000000, places=6)
        self.assertAlmostEqual(results[1]['DEITP_feature_deTP_mean_score'], 0.000000, places=6)
        self.assertAlmostEqual(results[1]['DEITP_feature_deTP_sd_score'], 0.000000, places=6)
        

    def test_float_deitp_yuv422p(self):
        ref_path, dis_path, asset, asset_original = set_default_576_324_yuv422p_videos_for_testing()
        self.fextractor = FloatDeitpFeatureExtractor(
            [asset, asset_original],
            None, fifo_mode=False,
            result_store=None
        )
        self.fextractor.run()
        results = self.fextractor.results

        #print(results[0])
        # print(results[1])

        self.assertAlmostEqual(results[0]['DEITP_feature_deitp_max_score'], 143.996747, places=5)
        self.assertAlmostEqual(results[0]['DEITP_feature_deitp_mean_score'], 14.587127, places=5)
        self.assertAlmostEqual(results[0]['DEITP_feature_deitp_sd_score'], 9.845331, places=5)
        self.assertAlmostEqual(results[1]['DEITP_feature_deitp_max_score'], 0.000000, places=6)
        self.assertAlmostEqual(results[1]['DEITP_feature_deitp_mean_score'], 0.000000, places=6)
        self.assertAlmostEqual(results[1]['DEITP_feature_deitp_sd_score'], 0.000000, places=6)

        self.assertAlmostEqual(results[0]['DEITP_feature_deTP_max_score'], 74.082346, places=5)
        self.assertAlmostEqual(results[0]['DEITP_feature_deTP_mean_score'], 5.932095, places=5)
        self.assertAlmostEqual(results[0]['DEITP_feature_deTP_sd_score'], 3.832568, places=5)
        self.assertAlmostEqual(results[1]['DEITP_feature_deTP_max_score'], 0.000000, places=6)
        self.assertAlmostEqual(results[1]['DEITP_feature_deTP_mean_score'], 0.000000, places=6)
        self.assertAlmostEqual(results[1]['DEITP_feature_deTP_sd_score'], 0.000000, places=6)


    def test_float_deitp_yuv422p10le(self):
        ref_path, dis_path, asset, asset_original = set_default_576_324_10bit_videos_for_testing()
        self.fextractor = FloatDeitpFeatureExtractor(
            [asset, asset_original],
            None, fifo_mode=False,
            result_store=None
        )
        self.fextractor.run()
        results = self.fextractor.results

        #print(results[0])
        #print(results[1])

        self.assertAlmostEqual(results[0]['DEITP_feature_deitp_max_score'], 145.094341, places=5)
        self.assertAlmostEqual(results[0]['DEITP_feature_deitp_mean_score'], 14.737297, places=5)
        self.assertAlmostEqual(results[0]['DEITP_feature_deitp_sd_score'], 9.866041, places=5)
        self.assertAlmostEqual(results[1]['DEITP_feature_deitp_max_score'], 0.000000, places=6)
        self.assertAlmostEqual(results[1]['DEITP_feature_deitp_mean_score'], 0.000000, places=6)
        self.assertAlmostEqual(results[1]['DEITP_feature_deitp_sd_score'], 0.000000, places=6)

        self.assertAlmostEqual(results[0]['DEITP_feature_deTP_max_score'], 79.513373, places=5)
        self.assertAlmostEqual(results[0]['DEITP_feature_deTP_mean_score'], 6.172106, places=5)
        self.assertAlmostEqual(results[0]['DEITP_feature_deTP_sd_score'], 4.026640, places=5)
        self.assertAlmostEqual(results[1]['DEITP_feature_deTP_max_score'], 0.000000, places=6)
        self.assertAlmostEqual(results[1]['DEITP_feature_deTP_mean_score'], 0.000000, places=6)
        self.assertAlmostEqual(results[1]['DEITP_feature_deTP_sd_score'], 0.000000, places=6)

        
    
    #Test Ipt input, BT2020 (assume yuv22p10le as ipt)
    def test_float_deitp_ipt422p10le_bt2020(self):
        ref_path, dis_path, asset, asset_original = set_default_576_324_10bit_videos_for_testing()
        self.fextractor = FloatDeitpFeatureExtractor(
            [asset, asset_original],
            None, fifo_mode=False,
            result_store=None,
            optional_dict={'eotf': 1, 'gamma': 2.4, 'clrFmt': 3, 'minLum': 0.005,
                           'maxLum': 100.0, 'Yuv2RgbExt': 'false',
                           'YuvXferSpec': 3, 'Rng': 0, 'RgbDef': 3}
            #optional_dict={'eotf': 0, 'gamma': 2.4}
        )
        self.fextractor.run()
        results = self.fextractor.results

        #print(results[0])
        #print(results[1])

        self.assertAlmostEqual(results[0]['DEITP_feature_deitp_max_score'], 286.888350, places=5)
        self.assertAlmostEqual(results[0]['DEITP_feature_deitp_mean_score'], 21.227075, places=5)
        self.assertAlmostEqual(results[0]['DEITP_feature_deitp_sd_score'], 15.620196, places=5)
        self.assertAlmostEqual(results[1]['DEITP_feature_deitp_max_score'], 0.000000, places=6)
        self.assertAlmostEqual(results[1]['DEITP_feature_deitp_mean_score'], 0.000000, places=6)
        self.assertAlmostEqual(results[1]['DEITP_feature_deitp_sd_score'], 0.000000, places=6)

        self.assertAlmostEqual(results[0]['DEITP_feature_deTP_max_score'], 93.238404, places=5)
        self.assertAlmostEqual(results[0]['DEITP_feature_deTP_mean_score'], 8.174942, places=5)
        self.assertAlmostEqual(results[0]['DEITP_feature_deTP_sd_score'], 5.496500, places=5)
        self.assertAlmostEqual(results[1]['DEITP_feature_deTP_max_score'], 0.000000, places=6)
        self.assertAlmostEqual(results[1]['DEITP_feature_deTP_mean_score'], 0.000000, places=6)
        self.assertAlmostEqual(results[1]['DEITP_feature_deTP_sd_score'], 0.000000, places=6)
    

    #Test Ipt input, BT2020 (assume yuv22p10le as ictcp)
    def test_float_deitp_ictcp422p10le_bt2020_fullrng(self):
        ref_path, dis_path, asset, asset_original = set_default_576_324_10bit_videos_for_testing()
        self.fextractor = FloatDeitpFeatureExtractor(
            [asset, asset_original],
            None, fifo_mode=False,
            result_store=None,
            optional_dict={'eotf': 1, 'gamma': 2.4, 'clrFmt': 4, 'minLum': 0.005,
                        'maxLum': 100.0, 'Yuv2RgbExt': 'false',
                        'YuvXferSpec': 3, 'Rng': 1, 'RgbDef': 3}
            #optional_dict={'eotf': 0, 'gamma': 2.4}
        )
        self.fextractor.run()
        results = self.fextractor.results

        #print(results[0])
        #print(results[1])

        self.assertAlmostEqual(results[0]['DEITP_feature_deitp_max_score'], 243.689144, places=5)
        self.assertAlmostEqual(results[0]['DEITP_feature_deitp_mean_score'], 17.739061, places=5)
        self.assertAlmostEqual(results[0]['DEITP_feature_deitp_sd_score'], 13.310461, places=5)
        self.assertAlmostEqual(results[1]['DEITP_feature_deitp_max_score'], 0.000000, places=6)
        self.assertAlmostEqual(results[1]['DEITP_feature_deitp_mean_score'], 0.000000, places=6)
        self.assertAlmostEqual(results[1]['DEITP_feature_deitp_sd_score'], 0.000000, places=6)

        self.assertAlmostEqual(results[0]['DEITP_feature_deTP_max_score'], 90.122321, places=5)
        self.assertAlmostEqual(results[0]['DEITP_feature_deTP_mean_score'], 6.173999, places=5)
        self.assertAlmostEqual(results[0]['DEITP_feature_deTP_sd_score'], 3.949051, places=5)
        self.assertAlmostEqual(results[1]['DEITP_feature_deTP_max_score'], 0.000000, places=6)
        self.assertAlmostEqual(results[1]['DEITP_feature_deTP_mean_score'], 0.000000, places=6)
        self.assertAlmostEqual(results[1]['DEITP_feature_deTP_sd_score'], 0.000000, places=6)


    def test_float_deitp_yuv444p(self):
        ref_path, dis_path, asset, asset_original = set_default_576_324_yuv444p_videos_for_testing()

        self.fextractor = FloatDeitpFeatureExtractor(
            [asset, asset_original],
            None, fifo_mode=False,
            result_store=None
        )
        self.fextractor.run()
        results = self.fextractor.results
        #print(results[0])
        # print(results[1])

        self.assertAlmostEqual(results[0]['DEITP_feature_deitp_max_score'], 142.677581, places=5)
        self.assertAlmostEqual(results[0]['DEITP_feature_deitp_mean_score'], 14.360697, places=5)
        self.assertAlmostEqual(results[0]['DEITP_feature_deitp_sd_score'], 9.803129, places=5)
        self.assertAlmostEqual(results[1]['DEITP_feature_deitp_max_score'], 0.000000, places=6)
        self.assertAlmostEqual(results[1]['DEITP_feature_deitp_mean_score'], 0.000000, places=6)
        self.assertAlmostEqual(results[1]['DEITP_feature_deitp_sd_score'], 0.000000, places=6)
        
        self.assertAlmostEqual(results[0]['DEITP_feature_deTP_max_score'], 67.212256, places=5)
        self.assertAlmostEqual(results[0]['DEITP_feature_deTP_mean_score'], 5.624041, places=5)
        self.assertAlmostEqual(results[0]['DEITP_feature_deTP_sd_score'], 3.576123, places=5)
        self.assertAlmostEqual(results[1]['DEITP_feature_deTP_max_score'], 0.000000, places=6)
        self.assertAlmostEqual(results[1]['DEITP_feature_deTP_mean_score'], 0.000000, places=6)
        self.assertAlmostEqual(results[1]['DEITP_feature_deTP_sd_score'], 0.000000, places=6)

    def test_float_deitp_rgb(self):
        ref_path, dis_path, asset, asset_original = set_default_576_324_yuv444p_videos_for_testing()

        # Check RBG color format (clrFmt = 1) consider same yuv input as rgb for time being
        self.fextractor = FloatDeitpFeatureExtractor(
            [asset, asset_original],
            None, fifo_mode=False,
            result_store=None,
            optional_dict={'eotf': 0, 'gamma': 2.4, 'clrFmt': 1, 'minLum': 0.005,
                           'maxLum': 100.0, 'Yuv2RgbExt': 'false',
                           'YuvXferSpec': 2, 'Rng': 0, 'RgbDef': 2}
            #optional_dict={'eotf': 0, 'gamma': 2.4}
        )

        self.fextractor.run()
        results = self.fextractor.results

        #print(results[0])
        #print(results[1])

        self.assertAlmostEqual(results[0]['DEITP_feature_deitp_max_score'], 90.046598, places=5)
        self.assertAlmostEqual(results[0]['DEITP_feature_deitp_mean_score'], 4.728545, places=5)
        self.assertAlmostEqual(results[0]['DEITP_feature_deitp_sd_score'], 4.166531, places=5)
        self.assertAlmostEqual(results[1]['DEITP_feature_deitp_max_score'], 0.000000, places=6)
        self.assertAlmostEqual(results[1]['DEITP_feature_deitp_mean_score'], 0.000000, places=6)
        self.assertAlmostEqual(results[1]['DEITP_feature_deitp_sd_score'], 0.000000, places=6)

        self.assertAlmostEqual(results[0]['DEITP_feature_deTP_max_score'], 85.345943, places=5)
        self.assertAlmostEqual(results[0]['DEITP_feature_deTP_mean_score'], 3.780756, places=5)
        self.assertAlmostEqual(results[0]['DEITP_feature_deTP_sd_score'], 3.983471, places=5)
        self.assertAlmostEqual(results[1]['DEITP_feature_deTP_max_score'], 0.000000, places=6)
        self.assertAlmostEqual(results[1]['DEITP_feature_deTP_mean_score'], 0.000000, places=6)
        self.assertAlmostEqual(results[1]['DEITP_feature_deTP_sd_score'], 0.000000, places=6)



    def test_float_deitp_yuv444p10le(self):
        ref_path, dis_path, asset, asset_original = set_default_576_324_yuv444p10le_videos_for_testing()
        self.fextractor = FloatDeitpFeatureExtractor(
            [asset, asset_original],
            None, fifo_mode=False,
            result_store=None
        )
        self.fextractor.run()
        results = self.fextractor.results

        #print(results[0])
        #print(results[1])

        self.assertAlmostEqual(results[0]['DEITP_feature_deitp_max_score'], 142.599217, places=5)
        self.assertAlmostEqual(results[0]['DEITP_feature_deitp_mean_score'], 14.307877, places=5)
        self.assertAlmostEqual(results[0]['DEITP_feature_deitp_sd_score'], 9.814519, places=5)
        self.assertAlmostEqual(results[1]['DEITP_feature_deitp_max_score'], 0.000000, places=6)
        self.assertAlmostEqual(results[1]['DEITP_feature_deitp_mean_score'], 0.000000, places=6)
        self.assertAlmostEqual(results[1]['DEITP_feature_deitp_sd_score'], 0.000000, places=6)

        self.assertAlmostEqual(results[0]['DEITP_feature_deTP_max_score'], 67.004993, places=5)
        self.assertAlmostEqual(results[0]['DEITP_feature_deTP_mean_score'], 5.529649, places=5)
        self.assertAlmostEqual(results[0]['DEITP_feature_deTP_sd_score'], 3.521925, places=5)
        self.assertAlmostEqual(results[1]['DEITP_feature_deTP_max_score'], 0.000000, places=6)
        self.assertAlmostEqual(results[1]['DEITP_feature_deTP_mean_score'], 0.000000, places=6)
        self.assertAlmostEqual(results[1]['DEITP_feature_deTP_sd_score'], 0.000000, places=6)


if __name__ == '__main__':
    unittest.main(verbosity=2)
