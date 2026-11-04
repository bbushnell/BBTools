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

/** V3 candidate-site vectors from measured mitochondrial raw anchors.
 * Training tables consume only TRAIN hosts. Every input row retains a disposition;
 * missing calls or unreachable ends never silently become correct examples.
 * This256-member cost-pilot tool does not modify the caller or fit a network.
 * Accepted-anchor rows are a subset of construction members; this tool is NOT
 * the all-member final-training pipeline. That pipeline must also align members
 * missing from, ambiguous in, or assigned to an original tier by the caller.
 * Training blocks hold the other raw end fixed; they are not sequential caller
 * predictions. Optional net/tables evaluation calls the actual sequential scorer
 * on original flanks. Neither mode installs a network in the production caller.
 * Preserved mitochondrial endpoint experiment: this class supports construction,
 * inference or validation of the experimental endpoint method. The method was not
 * promoted, and the mitochondrial-SSU direction was set aside on October 24 because
 * these genes did not fit the prokaryotic calling approach. Kept for a future
 * mitochondrial calling investigation; ordinary CallGenes never enables it.
 * @author Brian Bushnell, Raiden
 */
public final class MitoEndpointVectors {
	public static void main(String[] args){
		final PreParser pp=new PreParser(args,MitoEndpointVectors.class,false);
		try{final MitoEndpointVectors x=new MitoEndpointVectors();x.parse(pp.args);x.run();}finally{Shared.closeStream(pp.outstream);}
	}
	private void parse(String[] args){
		final Parser p=new Parser();
		for(String arg:args){final String[] kv=Parse.splitOnFirst(arg,'=');final String k=kv[0].toLowerCase(Locale.ROOT),v=kv.length>1?kv[1]:null;
			if(k.equals("anchors")){anchors=v;}else if(k.equals("outprefix")){prefix=v;}else if(k.equals("consensus")){consensus=v;}
			else if(k.equals("consensussha80")){consensusSha80=v;}else if(k.equals("expected")){expected=Parse.parseIntKMG(v);}
			else if(k.equals("startmin")){startMin=Integer.parseInt(v);}else if(k.equals("startmax")){startMax=Integer.parseInt(v);}
			else if(k.equals("stopmin")){stopMin=Integer.parseInt(v);}else if(k.equals("stopmax")){stopMax=Integer.parseInt(v);}
			else if(k.equals("startinside")){startInside=Integer.parseInt(v);}else if(k.equals("stopinside")){stopInside=Integer.parseInt(v);}
			else if(k.equals("net")){net=v;}else if(k.equals("tables")){savedTables=v;}
			else if(k.equals("ledger")){ledger=v;}
			else if(!p.parse(arg,k,v)){throw new IllegalArgumentException("Unknown flag: "+arg);}}
		fasta=p.in1;require(fasta!=null && anchors!=null && ledger!=null && prefix!=null && consensus!=null && SeedOffsetTable.validSha80(consensusSha80)
			&& expected==256 && !p.append && p.maxReads<0,"This cost-only256-member pilot cannot substitute for all-member final training");
		require(validRange(startMin,startMax) && validRange(stopMin,stopMax),"Explicit measured correction ranges must include0 and be bounded within +/-4096");
		Shared.TRIM_READ_DESCRIPTION=false;stream.Read.TO_UPPER_CASE=true;ReadWrite.USE_BGZF=true;ReadWrite.ALLOW_NATIVE_BGZF=true;ReadWrite.FORCE_BGZIP=true;ReadWrite.ZIPLEVEL=9;ReadWrite.setZipThreads(1);
		require((net==null)==(savedTables==null),"Evaluation must bind both saved network and training-only count tables");
		if(net==null){RrnaResourceIO.fresh(prefix,".train.tsv.gz",".validation.tsv.gz",".blocks.tsv.gz",".tables.tsv.gz");}
		else{RrnaResourceIO.fresh(prefix,".evaluation.tsv.gz");}
	}
	private void run(){
		final byte[][] refs=TrnaConsensusBuilder.loadLibrary(consensus);final String[] names=TrnaConsensusBuilder.lastLoadedNames;
		require(refs.length==158 && MitoSsuModelPolicy.coverage("16S",names),"Use the unchanged153+5 shipped consensus library");
		load(names);
		final float mean=MitoEndpointTables.medianLength(refs,5);
		if(net!=null){evaluate(refs,mean);return;}
		final MitoEndpointTables tables=new MitoEndpointTables(refs.length,startInside,stopInside,mean);
		for(Row r:rows){if(r.eligible && r.train){tables.add(r.bases,r.ts,r.te,r.model);}}
		tables.write(prefix+".tables.tsv.gz",consensusSha80);
		final ByteStreamWriter train=RrnaResourceIO.writer(prefix+".train.tsv.gz"),val=RrnaResourceIO.writer(prefix+".validation.tsv.gz"),blocks=RrnaResourceIO.writer(prefix+".blocks.tsv.gz");final ByteBuilder b=new ByteBuilder();
		final long[] vectors=new long[2];final int[] loci=new int[2];
		try{train.println("#dims\t10\t1");val.println("#dims\t10\t1");
			blocks.println("locus\thost\tsplit\tmodel_index\tstatus\tstart_rows\tstop_rows\ttruth_start0\ttruth_stop0\tanchor_start0\tanchor_stop0\tstart_min\tstart_max\tstop_min\tstop_max\tstart_inside\tstop_inside\tmean_length\tidentity\tgc");
			for(Row r:rows){
				int ns=0,ne=0;String status=r.status;
				if(r.eligible){
					final int c5=MitoEndpointGeometry.correction(r.ts,r.s),c3=MitoEndpointGeometry.correction(r.te,r.e);
					if(c5<startMin || c5>startMax || c3<stopMin || c3>stopMax){status="TRUTH_OUTSIDE_DECLARED_RANGE";}
					else if(r.e-r.ts<15 || r.te-r.s<15){status="OTHER_RAW_END_PREVENTS_TRUE_CANDIDATE";}
					else{
						final ByteStreamWriter out=r.train?train:val;
						ns=emit(r,tables,out,b,0,startMin,startMax);ne=emit(r,tables,out,b,1,stopMin,stopMax);
						status="VECTORS";final int split=r.train?0:1;vectors[split]+=ns+ne;loci[split]++;
					}
				}
				blocks.print(b.clear().append(r.id).tab().append(r.host).tab().append(r.train?"TRAIN":"VALIDATION").tab().append(r.model).tab().append(status)
					.tab().append(ns).tab().append(ne).tab().append(r.ts).tab().append(r.te).tab().append(r.s).tab().append(r.e)
					.tab().append(startMin).tab().append(startMax).tab().append(stopMin).tab().append(stopMax).tab().append(startInside).tab().append(stopInside)
					.tab().appendSlow(mean).tab().appendSlow(r.identity).tab().append(r.gc,8).nl());
			}
		}finally{RrnaResourceIO.finish(train,val,blocks);}
		require(loci[0]>0 && loci[1]>0,"Both host-separated splits need explicit eligible pilot loci before either trainer runs");
		System.err.println("MITO_ENDPOINT_VECTORS_PASS input_loci="+rows.size()+" train_loci="+loci[0]+" validation_loci="+loci[1]+" train_rows="+vectors[0]+" validation_rows="+vectors[1]+" output_files=4 k=9 features=v3 mean_length="+mean+" length_policy=mitochondrial_library_median independent_raw_blocks=true no_fit=true");
	}
	/** Harness check only: original raw anchors, confidence-ordered sequential
	 * sweeps, no gene discovery/path selection and no claim of net-on accuracy. */
	private void evaluate(byte[][] refs,float mean){
		final MitoEndpointTables tables=MitoEndpointTables.load(savedTables,refs.length,consensusSha80);
		require(tables.meanLength==mean && tables.insideCount(0)==startInside && tables.insideCount(1)==stopInside,
			"Evaluation must preserve the saved training feature definition");
		final CellNet model=NcrnaBoundaryScorer.load(net);require(model.numOutputs()==1,"Candidate-site confidence requires one output");
		final int[] starts=offsets(startMin,startMax),stops=offsets(stopMin,stopMax);
		final ByteStreamWriter out=RrnaResourceIO.writer(prefix+".evaluation.tsv.gz");final ByteBuilder b=new ByteBuilder();final int[] counts=new int[4];
		try{out.println("locus\thost\tsplit\tstatus\tmodel_index\tstart0\tstop0\texact5\texact3\tboth_exact\tstart_sites\tstop_sites\tclaim");
			for(Row r:rows){b.clear().append(r.id).tab().append(r.host).tab().append(r.train?"TRAIN":"VALIDATION").tab().append(r.status).tab().append(r.model);
				if(!r.eligible || r.train){out.print(b.append("\tNA\tNA\tNA\tNA\tNA\t0\t0\tNOT_EVALUATED").nl());continue;}
				final NcrnaBoundaryScorer.RefinementResult result=tables.refine(model,r.bases,r.s,r.e,r.model,refs[r.model],r.gc,r.identity,starts,stops,r.edgeMask);
				final int s=r.s+result.startOffset,e=r.e+result.stopOffset;final boolean a=s==r.ts,z=e==r.te;
				counts[0]++;if(a){counts[1]++;}if(z){counts[2]++;}if(a && z){counts[3]++;}
				out.print(b.tab().append(s).tab().append(e).tab().append(a).tab().append(z).tab().append(a && z)
					.tab().append(result.startSitesScored).tab().append(result.stopSitesScored).append("\tDEVELOPMENT_SCORER_HARNESS_NOT_CALLER_ACCURACY").nl());
			}
		}finally{RrnaResourceIO.finish(out);}
		require(counts[0]>0,"Sequential harness requires held-out-host anchors");
		System.err.println("MITO_ENDPOINT_HARNESS_PASS all_rows="+rows.size()+" validation_mito_anchors="+counts[0]+" exact5="+counts[1]+" exact3="+counts[2]+" both="+counts[3]+" caller_integration=false");
	}
	private static int[] offsets(int min,int max){
		require(validRange(min,max),"Scorer offsets must include the unchanged anchor");final int[] a=new int[max-min+1];for(int i=0;i<a.length;i++){a[i]=min+i;}return a;
	}
	private void load(String[] names){
		final ObjectMap<String,Row> map=new ObjectMap<String,Row>(String.class,Row.class);final ObjectMap<String,String> splits=new ObjectMap<String,String>(String.class,String.class);
		final ByteFile in=RrnaResourceIO.open(anchors);final LineParser1 p=new LineParser1('\t');
		try{RrnaResourceIO.header(in,ANCHOR_HEADER,anchors);
			for(byte[] line;(line=in.nextLine())!=null;){p.set(line);require(p.terms()==18,"Anchor ledger schema changed");final Row r=new Row();r.id=p.parseString(0);r.host=p.parseString(2);
				require(p.termEquals("TRAIN",3) || p.termEquals("VALIDATION",3),"Anchor split must remain a source-host split");r.train=p.termEquals("TRAIN",3);r.status=p.parseString(5);r.model=p.parseInt(8);r.ts=p.parseInt(10);r.te=p.parseInt(11);r.gc=p.parseFloat(17);
				require(r.model>=-1 && r.model<names.length,"Anchor model index must remain in the pinned library");
				final String split=p.parseString(3),old=splits.put(r.host,split);require(old==null || old.equals(split),"A host cannot contribute both training and validation loci");
				r.eligible=r.status.equals("MITO_EXACT_SOURCE_ANCHOR");
				if(r.model>=0){require(r.model<names.length && p.termEquals(names[r.model],9),"Anchor must bind the exact consensus name/index");r.s=p.parseInt(12);r.e=p.parseInt(13);r.identity=p.parseFloat(16);}
				if(r.eligible){require(p.termEquals("EXACT_SOURCE_ENDS",4) && p.parseInt(7)==1 && r.model>=5
					&& p.parseInt(14)==MitoEndpointGeometry.correction(r.ts,r.s) && p.parseInt(15)==MitoEndpointGeometry.correction(r.te,r.e),"Training eligibility must retain the measured raw mitochondrial anchor and correction sign");}
				require(map.put(r.id,r)==null,"One anchor ledger row per selected original member");rows.add(r);}
		}finally{require(!in.close(),"Anchor ledger read failed");}
		require(rows.size()==expected,"Never narrow the pilot denominator to only accepted training anchors");
		final ByteFile source=RrnaResourceIO.open(ledger);int sourceRows=0;
		try{RrnaResourceIO.header(source,MitoEndpointCorpus.HEADER,ledger);
			for(byte[] line;(line=source.nextLine())!=null;){p.set(line);require(p.terms()==18 && p.termEquals("true",12),"Window provenance must cover exactly the selected pilot");
				final Row r=map.get(p.parseString(0));require(r!=null && !r.provenance,"Every anchor needs exactly one source-window row");r.provenance=true;sourceRows++;
				require(p.termEquals(r.host,3) && p.termEquals(r.train?"TRAIN":"VALIDATION",4) && p.parseInt(14)==r.ts && p.parseInt(15)==r.te,
					"Source-window provenance must preserve host split and true endpoints");
				r.windowLength=p.parseInt(16);
				if(!p.termEquals("UNSUPPORTED_WINDOW",13)){
					require(p.termEquals("true",11) || p.termEquals("false",11),"Origin state must remain explicit");
					r.edgeMask=MitoEndpointTables.windowEdgeMask(new MitoEndpointGeometry(p.parseString(8),p.parseString(7),p.termEquals("true",11)),p.parseInt(6),r.ts,r.te,r.windowLength);
				}
			}
		}finally{require(!source.close(),"Source-window ledger read failed");}
		require(sourceRows==expected,"Source-window ledger cannot drop any selected member");
		RrnaResourceIO.stream(fasta,read->{final Row r=map.get(RrnaResourceIO.first(read.id));require(r!=null && r.bases==null,"Flanks must bind unique pilot rows");r.bases=read.bases;
			require(read.length()==r.windowLength && Integer.parseInt(RrnaResourceIO.field(read.id,"start0"))==r.ts && Integer.parseInt(RrnaResourceIO.field(read.id,"stop0"))==r.te
				&& RrnaResourceIO.field(read.id,"host").equals(r.host) && RrnaResourceIO.field(read.id,"split").equals(r.train?"TRAIN":"VALIDATION"),"Flanks must retain original truth and whole-host split");
			if(r.eligible){require(RrnaResourceIO.field(read.id,"partial").equals("false") && RrnaResourceIO.field(read.id,"pseudo").equals("false"),"Uncertain source ends cannot enter endpoint training");}});
		for(Row r:rows){if(r.eligible){require(r.bases!=null && r.ts>=0 && r.te<r.bases.length && r.s>=0 && r.e<r.bases.length
			&& r.s<=r.e && Float.isFinite(r.gc) && r.gc>=0 && r.gc<=1 && Float.isFinite(r.identity) && r.identity>=0 && r.identity<=1,"Eligible vectors require finite actual genomic context and accepted identity");}}
	}
	private static int emit(Row r,MitoEndpointTables tables,ByteStreamWriter out,ByteBuilder b,int end,int min,int max){
		assert(r.eligible && r.bases!=null):"Only explicit complete-source mito anchors supply endpoint labels";
		int rows=0,positives=0;
		for(int offset=min;offset<=max;offset++){
			final int s=r.s+(end==0?offset:0),e=r.e+(end==1?offset:0);if(s<0 || e>=r.bases.length || e-s<15){continue;}
			final int truth=end==0?r.ts:r.te;final float[] features=tables.features(r.bases,s,e,r.model,end,r.gc,r.identity,r.train?truth:-1);
			final int label=(end==0?s:e)==truth?1:0;positives+=label;b.clear();
			for(float f:features){require(Float.isFinite(f),"Nonfinite V3 inputs cannot reach a trainer");b.append(f,8).tab();}out.print(b.append(label).nl());rows++;
		}
		require(positives==1,"Every emitted endpoint block must retain exactly one reachable true endpoint");return rows;
	}
	private static boolean validRange(int min,int max){return min>=-4096 && min<=0 && max>=0 && max<=4096;}
	private static void require(boolean ok,String why){if(!ok){throw new IllegalArgumentException(why);}}
	private static final class Row{String id,host,status;byte[] bases;boolean train,eligible,provenance;int model,ts,te,s=-1,e=-1,edgeMask,windowLength;float gc,identity=Float.NaN;}
	private final ArrayList<Row> rows=new ArrayList<Row>();
	private String fasta,anchors,consensus,consensusSha80,prefix,net,savedTables,ledger;private int expected,startInside=-1,stopInside=-1,startMin=Integer.MIN_VALUE,startMax=Integer.MIN_VALUE,stopMin=Integer.MIN_VALUE,stopMax=Integer.MIN_VALUE;
	private static final String ANCHOR_HEADER="locus_id\tarm\thost\tsplit\tsource_status\tanchor_status\treturned_calls\toverlapping_calls\tmodel_index\tmodel\ttruth_start0\ttruth_stop0\tcaptured_start0\tcaptured_stop0\tcorrection5\tcorrection3\tidentity\tgc";
}
