package test.prok;

import prok.MitoSsuModelPolicy;
import prok.NcrnaBoundaryScorer;
import prok.RrnaResourceIO;
import prok.SeedOffsetTable;
import prok.TrnaConsensusBuilder;

import java.util.ArrayList;
import java.util.Locale;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.ReadWrite;
import map.ObjectMap;
import ml.CellNet;
import parse.LineParser1;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import shared.Shared;
import structures.ByteBuilder;

/** Full construction-ledger endpoint vectors and radius reachability.
 * Every member is retained in the locus ledger. Only complete associated raw
 * anchors supply training rows; unreachable true ends are reported explicitly.
 * Candidate blocks with zero reachable true endpoints remain labeled negatives,
 * not silently dropped members. No network is trained or installed here.
 * Preserved mitochondrial endpoint experiment: this class supports construction,
 * inference or validation of the experimental endpoint method. The method was not
 * promoted, and the mitochondrial-SSU direction was set aside on October 24 because
 * these genes did not fit the prokaryotic calling approach. Kept for a future
 * mitochondrial calling investigation; ordinary CallGenes never enables it.
 * @author Brian Bushnell, Raiden
 */
public final class MitoEndpointTrainingVectors {
	public static void main(String[] args){
		final PreParser pp=new PreParser(args,MitoEndpointTrainingVectors.class,false);
		try{final MitoEndpointTrainingVectors x=new MitoEndpointTrainingVectors();x.parse(pp.args);x.run();}finally{Shared.closeStream(pp.outstream);}
	}
	private void parse(String[] args){
		final Parser p=new Parser();
		for(String arg:args){final String[] kv=Parse.splitOnFirst(arg,'=');final String k=kv[0].toLowerCase(Locale.ROOT),v=kv.length>1?kv[1]:null;
			if(k.equals("assignments")){assignments=v;}else if(k.equals("ledger")){ledger=v;}else if(k.equals("consensus")){consensus=v;}
			else if(k.equals("consensussha80")){pin=v;}else if(k.equals("outprefix")){prefix=v;}
			else if(k.equals("startinside")){startInside=Integer.parseInt(v);}else if(k.equals("stopinside")){stopInside=Integer.parseInt(v);}
			else if(k.equals("net")){network=v;}else if(k.equals("tables")){savedTables=v;}else if(k.equals("radius")){radius=Integer.parseInt(v);radiusSet=true;}
			else if(k.equals("contextgenomes")){contextGenomes=v;}
			else if(k.equals("thresholds")){final String[] values=v.split(",");thresholds=new float[values.length];thresholdNames=new String[values.length];
				for(int i=0;i<values.length;i++){thresholds[i]=Float.parseFloat(values[i]);thresholdNames[i]=Float.toString(thresholds[i]);require(!Float.isNaN(thresholds[i]) && (i==0 || thresholds[i]>thresholds[i-1]),"Score thresholds must be strictly increasing and never NaN");}
				require(values.length>=2 && thresholds[0]==Float.NEGATIVE_INFINITY && thresholds[values.length-1]==Float.POSITIVE_INFINITY,"Gate sweeps require explicit ungated and all-raw controls");}
			else if(!p.parse(arg,k,v)){throw new IllegalArgumentException("Unknown flag: "+arg);}}
		fasta=p.in1;require(fasta!=null && assignments!=null && ledger!=null && consensus!=null && prefix!=null && SeedOffsetTable.validSha80(pin)
			&& !p.append && p.maxReads<0,"Full vector preparation requires pinned construction assignments, original genomic windows and the complete declared cohort");
		Shared.TRIM_READ_DESCRIPTION=false;stream.Read.TO_UPPER_CASE=true;
		ReadWrite.USE_BGZF=true;ReadWrite.ALLOW_NATIVE_BGZF=true;ReadWrite.FORCE_BGZIP=true;ReadWrite.ZIPLEVEL=9;ReadWrite.setZipThreads(1);
		require((network==null)==(savedTables==null) && (radius==128 || radius==192 || radius==256),
			"Evaluation binds both trained net and TRAIN-only tables at an explicit requested inference radius");
		require(network!=null || !radiusSet,"radius= controls evaluation only; full training-vector generation always retains both128/256datasets");
		require(network==null?contextGenomes==null:contextGenomes!=null,"Evaluation requires original contextgenomes= to expand artificial crops without changing raw anchors; vector construction remains unchanged");
		require(thresholds==null || network!=null,"Score gates apply only to saved-model evaluation, never training-vector construction");
		if(network==null){RrnaResourceIO.fresh(prefix,SUFFIXES);}
		else{RrnaResourceIO.fresh(prefix,".evaluation.tsv.gz",".evaluation.summary.tsv",".context.tsv");}
		if(thresholds!=null){RrnaResourceIO.fresh(prefix,".gates.tsv.gz",".confidence.tsv.gz");}
	}
	private void run(){
		final byte[][] refs=TrnaConsensusBuilder.loadLibrary(consensus);final String[] names=TrnaConsensusBuilder.lastLoadedNames;
		require(refs.length==158 && MitoSsuModelPolicy.coverage("16S",names),"Preserve the unchanged153mitochondrial plus five original models");
		load(names);
		if(network!=null){expandContexts();evaluate(refs);return;}
		final MitoEndpointTables tables=new MitoEndpointTables(refs.length,startInside,stopInside,MitoEndpointTables.medianLength(refs,5));
		for(Row r:rows){if(r.eligible && r.train){tables.add(r.bases,r.ts,r.te,r.model);}}
		tables.write(prefix+".tables.tsv.gz",pin);
		final ByteStreamWriter[][] writers=new ByteStreamWriter[2][2];
		for(int radius=0;radius<2;radius++){for(int split=0;split<2;split++){
			writers[radius][split]=RrnaResourceIO.writer(prefix+".r"+RADII[radius]+(split==0?".train.tsv.gz":".validation.tsv.gz"));writers[radius][split].println("#dims\t10\t1");
		}}
		final ByteStreamWriter detail=RrnaResourceIO.writer(prefix+".loci.tsv.gz");final ByteBuilder b=new ByteBuilder();
		final long[][][] counts=new long[CURVE_RADII.length][3][9];
		try{detail.println("locus\thost\tsplit\tconstruction_disposition\tvector_status\tmodel_index\twindow_edge_mask\tr128_start_rows\tr128_stop_rows\tr128_start_positives\tr128_stop_positives\tr256_start_rows\tr256_stop_rows\tr256_start_positives\tr256_stop_positives");
			for(Row r:rows){final int split=r.train?0:1;final long[][] c=new long[2][4],curve=new long[CURVE_RADII.length][4];String status=r.disposition;
				if(r.eligible){
					if(r.e-r.s<15){status="RAW_SPAN_TOO_SHORT_FOR_SCORER";}
					else{emit(r,tables,writers,b,c,curve);status="RAW_ANCHOR_CANDIDATE_BLOCKS";}
				}
				b.clear().append(r.id).tab().append(r.host).tab().append(r.train?"TRAIN":"VALIDATION").tab().append(r.disposition).tab().append(status).tab().append(r.model).tab().append(r.edgeMask);
				for(int radius=0;radius<2;radius++){for(long x:c[radius]){b.tab().append(x);}}
				for(int radius=0;radius<CURVE_RADII.length;radius++){accumulate(counts[radius][split],r,curve[radius]);accumulate(counts[radius][2],r,curve[radius]);}
				detail.print(b.nl());
			}
		}finally{RrnaResourceIO.finish(detail,writers[0][0],writers[0][1],writers[1][0],writers[1][1]);}
		final ByteStreamWriter summary=RrnaResourceIO.writer(prefix+".summary.tsv.gz");
		try{summary.println("radius\tsplit\tcohort_members\tbase_usable_members\tscored_members\tstart_reachable\tstop_reachable\tboth_reachable\tstart_candidate_rows\tstop_candidate_rows\tpositive_rows\tvector_bytes\tproposed_net_inputs\tproposed_net_parameters\tnet_weight_macs_per_candidate\tnet_activations_per_candidate\tmaximum_scores_per_locus\tclaim");
			for(int radius=0;radius<CURVE_RADII.length;radius++){for(int split=0;split<3;split++){
				final long[] c=counts[radius][split];require(c[0]==(split==0?12886:split==1?3224:16110) && c[1]==(split==0?12065:split==1?3019:15084),"Radius reporting must retain full-cohort and usable-anchor denominators");
				final int r=CURVE_RADII[radius];b.clear().append(r).tab().append(split==0?"TRAIN":split==1?"VALIDATION":"ALL");for(long x:c){b.tab().append(x);}b.tab();
				if(r==128 || r==256){b.append(vectorBytes(r==128?0:1,split));}else{b.append("NA_NOT_GENERATED");}
				summary.print(b.append("\t10\t193\t176\t17\t").append(4*r+4).append("\tREACHABILITY_AND_ANALYTIC_SCORING_WORK_NOT_NN_ACCURACY_OR_TIMED_LATENCY").nl());
			}}
		}finally{RrnaResourceIO.finish(summary);}
		System.err.println("MITO_ENDPOINT_TRAINING_VECTORS_PASS cohort=16110 usable=15084 train_usable=12065 validation_usable=3019 other_dispositions=1026 mean_length="+tables.meanLength+" radii=128,256 output_files=7 no_fit=true");
	}
	/** Score the exported network on every usable validation host, using the actual
	 * sequential endpoint scorer. This does not install it in CallGenes. */
	private void evaluate(byte[][] refs){
		final CellNet net=NcrnaBoundaryScorer.load(network);
		require(net!=null && net.numInputs()==10 && net.numOutputs()==1,"Shared sliding endpoint inference requires exactly10inputs/1output");
		final MitoEndpointTables tables=MitoEndpointTables.load(savedTables,refs.length,pin);
		require(tables.meanLength==MitoEndpointTables.medianLength(refs,5),"Evaluation must retain the trained family-wide length normalization");
		final int[] offsets=new int[2*radius+1];for(int i=0;i<offsets.length;i++){offsets[i]=i-radius;}
		final long[] totals=new long[9];final ByteBuilder b=new ByteBuilder();
		final ByteStreamWriter out=RrnaResourceIO.writer(prefix+".evaluation.tsv.gz");
		final ByteStreamWriter gates=thresholds==null?null:RrnaResourceIO.writer(prefix+".gates.tsv.gz"),confidence=thresholds==null?null:RrnaResourceIO.writer(prefix+".confidence.tsv.gz");
		long actualScores=0;
		try{out.println("locus\thost\tsplit\tdisposition\tmodel_index\tradius\traw_start0\traw_stop0\tstart0\tstop0\terror5\terror3\texact5\texact3\tboth_exact\tstart_sites\tstop_sites\tclaim");
			if(gates!=null){gates.println("locus\thost\tthreshold\terror5\terror3\texact5\texact3\tboth_exact\tmoved5\tmoved3\tchosen_score5\tchosen_score3\tclaim");
				confidence.println("locus\tstart_first\tinitial_start_score\tinitial_stop_score\tfirst_best\tfirst_raw\tsecond_moved_best\tsecond_moved_raw\tsecond_unmoved_best\tsecond_unmoved_raw\tphysical_network_scores\tfirst_best_offset\tsecond_moved_best_offset\tsecond_unmoved_best_offset\tfirst_best_float_bits\tfirst_raw_float_bits\tsecond_moved_best_float_bits\tsecond_moved_raw_float_bits\tsecond_unmoved_best_float_bits\tsecond_unmoved_raw_float_bits");}
			for(Row r:rows){b.clear().append(r.id).tab().append(r.host).tab().append(r.train?"TRAIN":"VALIDATION").tab().append(r.disposition).tab().append(r.model).tab().append(radius);
				if(!r.eligible || r.train){out.print(b.append("\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\t0\t0\tNOT_EVALUATED").nl());continue;}
				require(r.e-r.s>=15,"The reviewed full vector pass scored all3019usable validation anchors; evaluation cannot silently drop one");
				final NcrnaBoundaryScorer.RefinementResult result;
				MitoEndpointConfidenceGate cache=null;
				try{if(thresholds==null){result=tables.refine(net,r.bases,r.s,r.e,r.model,refs[r.model],r.gc,r.identity,offsets,offsets,r.edgeMask);}
					else{cache=new MitoEndpointConfidenceGate(net,r.bases,r.s,r.e,r.model,refs[r.model],r.gc,r.identity,tables,offsets,r.edgeMask);result=cache.at(Float.NEGATIVE_INFINITY);}}
				catch(AssertionError error){throw new IllegalArgumentException("Incomplete scoring for locus "+r.id+" host="+r.host+" raw="+r.s+".."+r.e+" truth="+r.ts+".."+r.te+" edges="+r.edgeMask,error);}
				require(Float.isFinite(result.startChosenScore) && Float.isFinite(result.stopChosenScore),"Nonfinite endpoint confidence is a failed model evaluation, not a correct unchanged endpoint");
				final int s=r.s+result.startOffset,e=r.e+result.stopOffset;
				final boolean a=s==r.ts,z=e==r.te,rawA=r.s==r.ts,rawZ=r.e==r.te;
				totals[0]++;if(rawA){totals[1]++;}if(rawZ){totals[2]++;}if(rawA && rawZ){totals[3]++;}
				if(a){totals[4]++;}if(z){totals[5]++;}if(a && z){totals[6]++;}totals[7]+=result.startSitesScored;totals[8]+=result.stopSitesScored;
				out.print(b.tab().append(r.s).tab().append(r.e).tab().append(s).tab().append(e).tab().append(s-r.ts).tab().append(e-r.te)
					.tab().append(a).tab().append(z).tab().append(a && z).tab().append(result.startSitesScored).tab().append(result.stopSitesScored)
					.append("\tVALIDATION_SEQUENTIAL_SCORER_NOT_WHOLE_GENOME_CALLER_ACCURACY").nl());
				if(cache!=null){actualScores+=cache.physicalScores;
					confidence.print(b.clear().append(r.id).tab().append(cache.startFirst).tab().append(cache.startConfidence,8).tab().append(cache.stopConfidence,8)
						.tab().append(cache.firstBest(),8).tab().append(cache.firstRaw(),8).tab().append(cache.secondMovedBest(),8).tab().append(cache.secondMovedRaw(),8)
						.tab().append(cache.secondUnmovedBest(),8).tab().append(cache.secondUnmovedRaw(),8).tab().append(cache.physicalScores)
						.tab().append(cache.firstBestOffset()).tab().append(cache.secondMovedBestOffset()).tab().append(cache.secondUnmovedBestOffset())
						.tab().append(Float.floatToIntBits(cache.firstBest())).tab().append(Float.floatToIntBits(cache.firstRaw()))
						.tab().append(Float.floatToIntBits(cache.secondMovedBest())).tab().append(Float.floatToIntBits(cache.secondMovedRaw()))
						.tab().append(Float.floatToIntBits(cache.secondUnmovedBest())).tab().append(Float.floatToIntBits(cache.secondUnmovedRaw())).nl());
					for(int i=0;i<thresholds.length;i++){final NcrnaBoundaryScorer.RefinementResult g=cache.at(thresholds[i]);final int d5=r.s+g.startOffset-r.ts,d3=r.e+g.stopOffset-r.te;
						gates.print(b.clear().append(r.id).tab().append(r.host).tab().append(thresholdNames[i]).tab().append(d5).tab().append(d3)
							.tab().append(d5==0).tab().append(d3==0).tab().append(d5==0 && d3==0).tab().append(g.startOffset!=0).tab().append(g.stopOffset!=0)
							.tab().append(g.startChosenScore,8).tab().append(g.stopChosenScore,8).append("\tSCORER_SEQUENTIAL_ABSOLUTE_SCORE_GATE_NOT_CALLER").nl());}
				}
			}
		}finally{RrnaResourceIO.finish(out);if(gates!=null){RrnaResourceIO.finish(gates,confidence);}}
		require(totals[0]==3019,"Validation model selection must account for every usable held-out-host anchor");
		final ByteStreamWriter summary=RrnaResourceIO.writer(prefix+".evaluation.summary.tsv");
		try{summary.println("radius\tvalidation_members\traw_exact5\traw_exact3\traw_both\tnet_exact5\tnet_exact3\tnet_both\tstart_scores\tstop_scores\tclaim");
			b.clear().append(radius);for(long n:totals){b.tab().append(n);}summary.print(b.append("\tVALIDATION_SEQUENTIAL_SCORER_NOT_WHOLE_GENOME_CALLER_ACCURACY").nl());
		}finally{RrnaResourceIO.finish(summary);}
		System.err.println("MITO_ENDPOINT_FULL_EVALUATION_PASS validation=3019 radius="+radius+" raw_both="+totals[3]+" net_both="+totals[6]+" caller_integration=false");
		if(thresholds!=null){System.err.println("MITO_ENDPOINT_CONFIDENCE_GATES_PASS validation=3019 thresholds="+thresholds.length+" rows="+(3019L*thresholds.length)+" physical_network_scores="+actualScores);}
	}
	private long vectorBytes(int radius,int split){
		require(radius>=0 && radius<2 && split>=0 && split<3,"Vector sizes must identify a declared radius and host split");
		if(split==2){return vectorBytes(radius,0)+vectorBytes(radius,1);}
		final java.io.File file=new java.io.File(prefix+".r"+RADII[radius]+(split==0?".train.tsv.gz":".validation.tsv.gz"));
		require(file.isFile() && file.length()>0,"Closed vector files must exist before their byte sizes enter a result");return file.length();
	}
	/** Extend the saved crop using its original genome, preserving every anchor,
	 * identity, host split and correction. No alignment or training is repeated. */
	private void expandContexts(){
		final ObjectMap<String,Rows> byAccession=new ObjectMap<String,Rows>(String.class,Rows.class);
		for(Row r:rows){if(r.eligible && !r.train){Rows list=byAccession.get(r.accession);if(list==null){list=new Rows();byAccession.put(r.accession,list);}list.add(r);}}
		final ByteStreamWriter out=RrnaResourceIO.writer(prefix+".context.tsv");final ByteBuilder b=new ByteBuilder();final int[] count={0};
		try{out.println("locus\taccession\toriginal_window_length\texpanded_window_length\tcoordinate_shift\tphysical_edge_mask\texplicit_origin_crossing");
			RrnaResourceIO.stream(contextGenomes,read->{final ArrayList<Row> list=byAccession.get(RrnaResourceIO.first(read.id));if(list==null){return;}
				for(Row r:list){require(!r.expanded && read.length()==r.genomeLength,"Original source identity/length must match each validation window once: "+r.id);
					final MitoEndpointGeometry.Window original=r.geometry.window(read.bases,1024,32768);
					require(original!=null && original.start==r.ts && original.stop==r.te && java.util.Arrays.equals(original.bases,r.bases),"Reconstructed original window must exactly equal the frozen corpus before adding context: "+r.id);
					// enrichmentProfileCore reads candidate +/-2. Across every legal
					// inside/outside placement, kmerWindowStart can extend k bases
					// toward either side; k+2 covers that full feature footprint.
					final MitoEndpointGeometry.Window expanded=r.geometry.extendWindow(read.bases,original,radius+MitoEndpointTables.K+2);
					final int shift=expanded.start-r.ts,oldLength=r.bases.length;
					r.ts+=shift;r.te+=shift;r.s+=shift;r.e+=shift;r.bases=expanded.bases;r.length=r.bases.length;r.expanded=true;count[0]++;
					r.edgeMask=r.geometry.crossesOrigin?0:MitoEndpointTables.windowEdgeMask(r.geometry,r.genomeLength,r.ts,r.te,r.length);
					out.print(b.clear().append(r.id).tab().append(r.accession).tab().append(oldLength).tab().append(r.length).tab().append(shift).tab().append(r.edgeMask).tab().append(r.geometry.crossesOrigin).nl());
				}
			});
		}finally{RrnaResourceIO.finish(out);}
		require(count[0]==3019,"Every held-out anchor requires checked original-genome context; none may be silently dropped");
	}
	private static void accumulate(long[] total,Row r,long[] c){
		assert(total.length==9 && c.length==4):"Report cohort, reachability and actual candidate rows separately";
		total[0]++;if(r.eligible){total[1]++;}if(c[0]>0 && c[1]>0){total[2]++;}
		total[3]+=c[2];total[4]+=c[3];if(c[2]==1 && c[3]==1){total[5]++;}total[6]+=c[0];total[7]+=c[1];total[8]+=c[2]+c[3];
	}
	private static void emit(Row r,MitoEndpointTables tables,ByteStreamWriter[][] out,ByteBuilder b,long[][] counts,long[][] curve){
		assert(r.eligible && r.bases!=null && r.e-r.s>=15):"Training rows require a supported raw candidate span, not merely an assigned model";
		writeCandidates(r.bases,r.s,r.e,r.ts,r.te,r.model,r.gc,r.identity,r.train,tables,out[0][r.train?0:1],out[1][r.train?0:1],b,counts,curve);
	}
	/** Shared row writer; a zero-positive block is retained and counted explicitly. */
	public static void writeCandidates(byte[] bases,int rawStart,int rawStop,int trueStart,int trueStop,int model,float gc,float identity,
			boolean train,MitoEndpointTables tables,ByteStreamWriter narrow,ByteStreamWriter wide,ByteBuilder b,long[][] counts){
		writeCandidates(bases,rawStart,rawStop,trueStart,trueStop,model,gc,identity,train,tables,narrow,wide,b,counts,null);
	}
	public static void writeCandidates(byte[] bases,int rawStart,int rawStop,int trueStart,int trueStop,int model,float gc,float identity,
			boolean train,MitoEndpointTables tables,ByteStreamWriter narrow,ByteStreamWriter wide,ByteBuilder b,long[][] counts,long[][] curve){
		require(bases!=null && rawStart>=0 && rawStop<bases.length && rawStop-rawStart>=15 && trueStart>=0 && trueStop<bases.length && trueStop>=trueStart
			&& counts.length==2 && counts[0].length==4 && counts[1].length==4,"Candidate blocks require valid original/raw coordinates and separate radius counters");
		for(int end=0;end<2;end++){
			final int truth=end==0?trueStart:trueStop;
			for(int offset=-256;offset<=256;offset++){
				final int s=rawStart+(end==0?offset:0),e=rawStop+(end==1?offset:0);if(s<0 || e>=bases.length || e-s<15){continue;}
				final float[] features=tables.features(bases,s,e,model,end,gc,identity,train?truth:-1);
				final int label=(end==0?s:e)==truth?1:0;b.clear();for(float value:features){require(Float.isFinite(value),"Every training feature must be finite");b.append(value,8).tab();}b.append(label).nl();
				wide.print(b);counts[1][end]++;counts[1][end+2]+=label;
				if(Math.abs(offset)<=128){narrow.print(b);counts[0][end]++;counts[0][end+2]+=label;}
				if(curve!=null){for(int i=0;i<CURVE_RADII.length;i++){if(Math.abs(offset)<=CURVE_RADII[i]){curve[i][end]++;curve[i][end+2]+=label;}}}
			}
			for(int i=0;i<2;i++){final int radius=RADII[i];final boolean reachable=Math.abs(truth-(end==0?rawStart:rawStop))<=radius
				&& (end==0?rawStop-trueStart:trueStop-rawStart)>=15;
				require(counts[i][end]>0 && counts[i][end+2]==(reachable?1:0),"Unreachable endpoints remain explicit zero-positive blocks; reachable ones require exactly one positive row");}
		}
	}
	private void load(String[] names){
		final ObjectMap<String,Row> map=new ObjectMap<String,Row>(String.class,Row.class);final ObjectMap<String,String> splits=new ObjectMap<String,String>(String.class,String.class);
		final ByteFile in=RrnaResourceIO.open(assignments);final LineParser1 p=new LineParser1('\t');int eligible=0,train=0,val=0;
		try{RrnaResourceIO.header(in,HEADER,assignments);
			for(byte[] line;(line=in.nextLine())!=null;){p.set(line);require(p.terms()==21,"Full construction-assignment schema changed");final Row r=new Row();r.id=p.parseString(0);r.host=p.parseString(2);
				require(p.termEquals("TRAIN",3) || p.termEquals("VALIDATION",3),"Host split must remain explicit");r.train=p.termEquals("TRAIN",3);
				final String split=p.parseString(3),old=splits.put(r.host,split);require(old==null || old.equals(split),"One host cannot supply both TRAIN and VALIDATION");
				r.disposition=p.parseString(6);r.model=p.parseInt(7);r.ts=p.parseInt(13);r.te=p.parseInt(14);r.s=p.parseInt(15);r.e=p.parseInt(16);
				require(r.model>=5 && r.model<names.length && p.termEquals(names[r.model],8),"Every full-cohort assignment must bind a mitochondrial model in the pinned library");
				r.eligible=r.disposition.equals("COMPLETE_SOURCE_ASSOCIATED");
				if(r.eligible){require(p.termEquals("EXACT_SOURCE_ENDS",4) && p.termEquals("false",19) && p.termEquals("false",20)
					&& p.parseInt(17)==r.ts-r.s && p.parseInt(18)==r.te-r.e,"Usable anchors require exact source ends, conserved correction signs and no uncertain-source flags");
					r.identity=p.parseFloat(12);require(Float.isFinite(r.identity) && r.identity>=0 && r.identity<=1,"V3 uses the recorded full-model genomic identity");eligible++;if(r.train){train++;}else{val++;}}
				require(map.put(r.id,r)==null,"Every construction member must occur exactly once");rows.add(r);}
		}finally{require(!in.close(),"Construction ledger read failed");}
		require(rows.size()==16110 && eligible==15084 && train==12065 && val==3019,"Full training preparation must preserve the reviewed cohort and usable-anchor split");
		final ByteFile source=RrnaResourceIO.open(ledger);int n=0;
		try{RrnaResourceIO.header(source,MitoEndpointCorpus.HEADER,ledger);
			for(byte[] line;(line=source.nextLine())!=null;){p.set(line);require(p.terms()==18 && p.termEquals("true",12),"Use every selected full-cohort source window");final Row r=map.get(p.parseString(0));
				require(r!=null && !r.provenance && p.termEquals(r.host,3) && p.termEquals(r.train?"TRAIN":"VALIDATION",4)
					&& p.parseInt(14)==r.ts && p.parseInt(15)==r.te,"Window provenance must match each original host split and true endpoint");r.provenance=true;n++;
				r.length=p.parseInt(16);r.genomeLength=p.parseInt(6);r.accession=p.parseString(2);r.gc=p.parseFloat(17);if(r.eligible){
					require(p.termEquals("false",9) && p.termEquals("false",10) && p.termEquals("EXACT_SOURCE_ENDS",13),"Uncertain source metadata cannot enter endpoint training");
					require(p.termEquals("true",11) || p.termEquals("false",11),"Origin crossing must remain an explicit source boolean");
					r.geometry=new MitoEndpointGeometry(p.parseString(8),p.parseString(7),p.termEquals("true",11));
					r.edgeMask=MitoEndpointTables.windowEdgeMask(r.geometry,r.genomeLength,r.ts,r.te,r.length);
				}
			}
		}finally{require(!source.close(),"Source-window ledger read failed");}
		require(n==16110,"Every construction row requires its original window disposition");
		final int[] windows={0};RrnaResourceIO.stream(fasta,read->{final Row r=map.get(RrnaResourceIO.first(read.id));require(r!=null && !r.windowSeen && read.length()==r.length,"Each original genomic window must occur exactly once");r.windowSeen=true;windows[0]++;
			require(Integer.parseInt(RrnaResourceIO.field(read.id,"start0"))==r.ts && Integer.parseInt(RrnaResourceIO.field(read.id,"stop0"))==r.te
				&& RrnaResourceIO.field(read.id,"host").equals(r.host) && RrnaResourceIO.field(read.id,"split").equals(r.train?"TRAIN":"VALIDATION"),"Actual flanks must retain the sealed source coordinates and split");if(r.eligible){r.bases=read.bases;}});
		require(windows[0]==16108,"Preserve both explicitly unsupported windows in the full-cohort accounting");
		for(Row r:rows){if(r.eligible){require(r.bases!=null && r.ts>=0 && r.te>=r.ts && r.te<r.length && r.s>=0 && r.e>=r.s && r.e<r.length && Float.isFinite(r.gc) && r.gc>=0 && r.gc<=1,"Usable vectors need actual finite genomic context and raw coordinates");}}
	}
	private static void require(boolean ok,String why){if(!ok){throw new IllegalArgumentException(why);}}
	private static final class Row{String id,host,disposition,accession;byte[] bases;int model,ts,te,s,e,length,edgeMask,genomeLength;float gc,identity;MitoEndpointGeometry geometry;boolean train,eligible,provenance,windowSeen,expanded;}
	private static final class Rows extends ArrayList<Row>{private static final long serialVersionUID=1L;}
	private final ArrayList<Row> rows=new ArrayList<Row>();private String fasta,assignments,ledger,consensus,pin,prefix,network,savedTables,contextGenomes;private int startInside=-1,stopInside=-1,radius=192;private boolean radiusSet;
	private float[] thresholds;private String[] thresholdNames;
	private static final int[] RADII={128,256};
	private static final int[] CURVE_RADII={32,48,64,96,128,192,256};
	private static final String[] SUFFIXES={".r128.train.tsv.gz",".r128.validation.tsv.gz",".r256.train.tsv.gz",".r256.validation.tsv.gz",".tables.tsv.gz",".loci.tsv.gz",".summary.tsv.gz"};
	private static final String HEADER="locus_id\tarm\thost\tsplit\tsource_status\tpilot_caller_status\tdisposition\tmodel_index\tmodel\tshared9\tmodels_aligned\tassignment_identity_shorter_query\tgenomic_identity_full_model\ttruth_start0\ttruth_stop0\traw_start0\traw_stop0\tcorrection5\tcorrection3\tpartial\tpseudo";
}
