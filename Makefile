INCDIR = -I.
DBG    = -g 
OPT    = -O3
CPP    = g++ 
CFLAGS = $(DBG) $(OPT) $(INCDIR) 
LINK   = -lm 

.cpp.o:
	$(CPP) $(CFLAGS) -c $< -o $@

all: segment # collate_sequences size_converter

segment: scene-labeling.cpp segment-image.h segment-video.h segment-graph.h disjoint-set.h label-components.h semantic-segment.h graph-matching.h list_utils.h color_flow.h basic-morpho.h
	$(CPP) $(CFLAGS) -o segment scene-labeling.cpp $(LINK)

size_converter: size_converter.cpp segment-image.h segment-video.h segment-graph.h disjoint-set.h label-components.h semantic-segment.h graph-matching.h list_utils.h color_flow.h basic-morpho.h
	$(CPP) $(CFLAGS) -o size_converter size_converter.cpp $(LINK)

collate_sequences: collate_sequences.cpp segment-image.h segment-video.h segment-graph.h disjoint-set.h label-components.h semantic-segment.h graph-matching.h list_utils.h color_flow.h basic-morpho.h
	$(CPP) $(CFLAGS) -o collate_sequences collate_sequences.cpp $(LINK)


clean:
	/bin/rm -f segment *.o *~

clean-all: clean
	/bin/rm  -f *~ out*.png; rm -f opticalout*.png; rm -f labeledout*.png; rm -f alone_segm*.png ; rm out*.ppm; rm labeledout*.ppm;


