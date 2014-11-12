import numpy as np
import h5py
import vigra
import pgmlink
import datetime

def getAnnotations(fn, path='/ManualTransitions/Labels/0', ):
   segment2Label = {}
   label2Segment = {}

   with h5py.File(fn, 'r') as f:
      labels = f[path]
      for t in labels.keys():
         labels_at = labels[t]
         t = int(t)
         if t not in segment2Label.keys():
            segment2Label[t] = {}
         if t not in label2Segment.keys():
            label2Segment[t] = {}

         for segment in labels_at.keys():
            segment2Label[t][int(float(segment))] = []
            for label in np.array(labels_at[segment]):
               segment2Label[t][int(float(segment))].append(int(float(label)))
               if int(float(label)) not in label2Segment[t].keys():
                  label2Segment[t][int(float(label))] = []
               label2Segment[t][int(float(label))].append(int(float(segment)))

   return segment2Label, label2Segment


def trainGPClassifier(featureArrays, labels, out_fn, out_path, feature_names):
    from lazyflow.classifiers.gaussianProcessClassifier import GaussianProcessClassifierFactory, GaussianProcessClassifier

    gcf = GaussianProcessClassifierFactory()
    gpc = gcf.create_and_train(featureArrays, labels)

    print 'writing to file'
    with h5py.File(out_fn, 'r+') as out:
        try:
            del out[out_path]
        except:
            pass
        g = out.create_group(out_path)
        gpc.serialize_hdf5(g)
        h = g.create_group('Features')
        for op, v in feature_names.items():
            h.create_dataset(data=v, name=op)
        g.create_dataset(name='timestamp', data=str(datetime.date.today()))

    return gpc


def computeFeatureStore(traxel_ids, labelImage, raw):
    print '  computing feature store...'
    ts = pgmlink.TraxelStore()
    fs = pgmlink.FeatureStore()

    for t in sorted(traxel_ids.keys()):
        for idx in traxel_ids[t]:
            trax = pgmlink.Traxel()
            trax.set_feature_store(fs)
            trax.Id = int(idx)
            trax.Timestep = int(t)
            ts.add(fs, trax)

    # compute features from labelimage / raw data
    for t in sorted(traxel_ids.keys()):
        li_at = labelImage[t,...,0].astype(np.uint32).squeeze()
        # TODO: relabel with traxel_ids???
        raw_at = raw[t,...,0].astype(np.float32).squeeze()
        print '    extracting features for timestep', t
        print '    ', raw_at.shape, li_at.shape, raw_at.dtype, li_at.dtype
        pgmlink.extract_region_features(raw_at, li_at, fs, t)

    return ts, fs


def initializeFeatureExtractors(featureNames):
    print '  initializing feature extractors...'
    extractors = []
    for operator in featureNames.keys():
        for featname in featureNames[operator]:
            extractors.append(pgmlink.FeatureExtractor(operator, featname))
    return extractors


def getTransitionClassifier(raw_fn, raw_path, fn, annotationPath, transitionFeatureNames, labelImg_fn,
                            labelImg_path, t_axis, ch_axis, ch, out_fn, out_path, margin=None,
                            max_distance=1000):
   print 'Transition Classifier:'
   print '======================'
   print '  reading annotations...'
   segment2Label, label2Segment = getAnnotations(fn, annotationPath)

   featureArrays = []
   # same ordering as featureArrays
   # label 1 means negative example, label 2 means positive example
   labels = []

   traxel_ids = {}
   for t in label2Segment.keys():
       for label in label2Segment[t].keys():
           #assert len(label2Segment[t][label]) == 1, \
           #      'only single segment examples are allowed; see the joint-tracking trainClassifier script for more'
           traxel_ids.setdefault(t, []).append(label2Segment[t][label][0])

   labelImage = h5py.File(labelImg_fn)[labelImg_path]
   raw = h5py.File(raw_fn)[raw_path]

   ts, fs = computeFeatureStore(traxel_ids, labelImage, raw)

   extractors = initializeFeatureExtractors(transitionFeatureNames)

   for t in label2Segment.keys():
      print '  processing timestep', t
      if t+1 not in label2Segment.keys():
         continue

      positivePairs_at = []
      negativePairs_at = []
      
      for label in label2Segment[t]:
         # if this label is not present in the next time step, this pair has already been added or is invalid

         if label not in label2Segment[t+1].keys():
            continue
         #assert len(label2Segment[t][label]) == 1, \
         #    'only single segment examples are allowed; see the joint-tracking trainClassifier script for more'
         pair = [ label2Segment[t][label][0], label2Segment[t+1][label][0] ]
         if label < 0:
            negativePairs_at.append(pair)
         elif label > 0:
            positivePairs_at.append(pair)
         else:
            raiseException, "labels must not be zero"
      
      if len(positivePairs_at) + len(negativePairs_at) == 0:
         # nothing to do in this time step
         continue

      for idx,pair in enumerate(positivePairs_at + negativePairs_at):
         trax1 = ts.get_traxel(pair[0], t)
         trax2 = ts.get_traxel(pair[1], t+1)

         # fs.print_traxel_features(t, pair[0])

         feats_at_idx = []
         assert transitionFeatures.keys() == ['SquaredDifference'] and \
             transitionFeatures['SquaredDifference'] == ['RegionCenter'], "max_distance only defined for SquaredDistance"
         for fe in extractors:
             v = fe.extract(trax1, trax2)
             if np.all(np.array(v) < max_distance):
                feats_at_idx.append(v)

         if len(feats_at_idx) == 0:
             continue

         feats_at_idx = np.array(feats_at_idx)

         if idx < len(positivePairs_at):
            labels.append(2)
         else:
            labels.append(1)

         featureArrays.append(feats_at_idx)

   featureArrays = np.array(featureArrays).squeeze()
   labels = np.array(labels)
   print '  training Gaussian process classifier...'
   gpc = trainGPClassifier(featureArrays, labels, out_fn, '/TransitionGPClassifier', transitionFeatureNames)
   print '  =================================='
   print '  Transition Classifier: done'
   print '  =================================='

   return gpc


if __name__ == '__main__':
   import sys
   argv = sys.argv
   name = 'mitocheck'
   if len(argv) > 1:
      name = argv[1]
   if name == 'mitocheck':
      fn = '/home/mschiegg/extern/data/cvpr2015/mitocheck/training/training_transition-classifier.ilp'
      raw_fn = '/home/mschiegg/extern/data/cvpr2015/mitocheck/training/training_raw.h5'
      pathRaw = 'exported_data'
      out_fn = '/home/mschiegg/extern/data/cvpr2015/mitocheck/training/training_transition-classifier.ilp'
      pathTransitions = '/ManualTransitions/Labels/0'

      labelImg_fn = '/home/mschiegg/extern/data/cvpr2015/mitocheck/training/training_labels.h5'
      pathLabelImg = 'exported_data_T'
      ch_axis = 3
      ch = 0
      t_axis = 0

   elif name == 'drosophila':
      data_dir = '/home/phanslov/ma/data/drosophila'
      fn = '%s/classifier_training.ilp' % data_dir
      raw_fn = '%s/subset_raw.h5' % data_dir
      pathRaw = 'volume/data'
      out_fn = '%s/drosophila_classifiers.h5' % data_dir
      pathTransitions = '/ManualTransitions/Labels/0'
      pathRegions = '/ManualRegions/Labels/0'
      pathDivisions = '/ManualDivisions/Labels/0'

      segmentImg_fn =  '%s/merged_oversegmentation.h5' % data_dir
      ch_axis = -1
      ch = 0
      segmentImg_path = 'exported_data'
      t_axis = 0

   elif name == 'drosophila_from_conservation':
      data_dir = '/home/mschiegg/extern/tmp'
      fn = '%s/conservationTracking_2013-08-23_transformed.ilp' % data_dir
      raw_fn = '%s/300-399_cropped_z30.h5' % data_dir
      pathRaw = 'volume/data'
      out_fn = '%s/drosophila_classifiers_division.h5' % data_dir
      pathTransitions = '/ManualTransitions/Labels/0'
      pathRegions = '/ManualRegions/Labels/0'
      pathDivisions = '/ManualDivisions/Labels/0'

      segmentImg_fn =  '%s/conservationTracking_labelImage.h5' % data_dir
      ch_axis = -1
      ch = 0
      segmentImg_path = 'exported_data'
      t_axis = 0
      div_from_ilastik=True

   elif name == 'rapoport1':
      fn = '/home/mschiegg/extern/data/cvpr2015/rapoport/dataset1/training1_transition-classifier.ilp'
      raw_fn = '/home/mschiegg/extern/data/cvpr2015/rapoport/dataset1/training1_raw.h5'
      pathRaw = 'exported_data'
      out_fn = '/home/mschiegg/extern/data/cvpr2015/rapoport/dataset1/training1_transition-classifier.ilp'
      pathTransitions = '/ManualTransitions/Labels/0'

      labelImg_fn = '/home/mschiegg/extern/data/cvpr2015/rapoport/dataset1/training1_labels.h5'
      pathLabelImg = 'exported_data_T'
      ch_axis = -1
      ch = 0
      segmentImg_path = 'exported_data'
      t_axis = 0
      
   elif name == 'rapoport2':
      fn = '/home/mschiegg/extern/data/cvpr2015/rapoport/dataset2/training2_transition-classifier.ilp'
      raw_fn = '/home/mschiegg/extern/data/cvpr2015/rapoport/dataset2/training2_raw.h5'
      pathRaw = 'exported_data_T'
      out_fn = '/home/mschiegg/extern/data/cvpr2015/rapoport/dataset2/training2_transition-classifier.ilp'
      pathTransitions = '/ManualTransitions/Labels/0'

      labelImg_fn = '/home/mschiegg/extern/data/cvpr2015/rapoport/dataset2/training2_labels.h5'
      pathLabelImg = 'exported_data'
      ch_axis = -1
      ch = 0
      segmentImg_path = 'exported_data'
      t_axis = 0
      
   else:
      raise Exception

   transitionFeatures = { 'SquaredDifference': ['RegionCenter'], }

   getTransitionClassifier(raw_fn, pathRaw, fn, pathTransitions, transitionFeatures, labelImg_fn, pathLabelImg, t_axis, ch_axis, ch, out_fn, 'Transitions')
