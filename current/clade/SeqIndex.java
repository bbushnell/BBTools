package clade;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;

import aligner.IDAligner;
import aligner.QuantumAligner;
import bin.BinObject;
import bin.BinSketcher;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import shared.Shared;
import sketch.CompareBuffer;
import sketch.Comparison;
import sketch.DisplayParams;
import sketch.Sketch;
import sketch.SketchIndex;
import sketch.SketchMakerMini;
import sketch.SketchObject;
import sketch.SketchResults;
import sketch.SketchTool;
import stream.ConcurrentReadInputStream;
import stream.Read;
import structures.ListNum;
import tax.TaxNode;
import tax.TaxTree;

public class SeqIndex extends CladeObject {
	
	public static void main(String[] args) {
		
//		String //TODO
		
		loadTree();
		
	}

	/*--------------------------------------------------------------*/
	/*----------------         Constructor          ----------------*/
	/*--------------------------------------------------------------*/

	public SeqIndex() {}

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/
	
	public long addFile(String fname, String type) {
		FileFormat ff=FileFormat.testInput(fname, FileFormat.FA, null, true, false);
		return addFile(ff, type);
	}
	
	public long addFile(FileFormat ff, String type) {
		ConcurrentReadInputStream cris=ConcurrentReadInputStream.getReadInputStream(-1, false, ff, null);
		cris.start();
		long added=add(cris, type);
		boolean errorState=ReadWrite.closeStreams(cris);
		if(verbose){System.err.println("Finished reading data.");}
		if(errorState){System.err.println("Something went wrong reading "+ff.name());}
		return added;
	}
	
	public long add(Read r, String type) {
		int tid=bin.BinObject.parseTaxID(r.name());
		if(tid<1) {return 0;}
		Organism o=index.get(tid);
		if(o==null) {
			TaxNode tn=tree.getNode(tid);
			if(tn==null) {return 0;}
			o=new Organism(tn.name, tid);
			index.put(tid, o);
			organismsAdded++;
		}
		if(maxSeqsPerOrganism>0 && o.seqs.size()>=maxSeqsPerOrganism) {return 0;}
		Sequence s=new Sequence(r.bases, r.name(), type, tid, r.numericID);
		o.add(s);
		sequencesAdded++;
		basesAdded+=s.size();
		return 1;
	}
	
	public long add(ConcurrentReadInputStream cris, String type) {
		ListNum<Read> ln=cris.nextList();
		ArrayList<Read> reads=(ln!=null ? ln.list : null);
		long added=0;
		while(ln!=null && reads!=null && reads.size()>0){
			for(int idx=0; idx<reads.size(); idx++){
				final Read r1=reads.get(idx), r2=r1.mate;
				readsProcessed+=r1.pairCount();
				basesProcessed+=r1.pairLength();
				added+=add(r1, type);
				if(r2!=null) {added+=add(r2, type);}
			}

			cris.returnList(ln);
			if(verbose){System.err.println("Returned a list.");}
			ln=cris.nextList();
			reads=(ln!=null ? ln.list : null);
		}
		if(ln!=null){
			cris.returnList(ln.id, ln.list==null || ln.list.isEmpty());
		}
		return added;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fetch             ----------------*/
	/*--------------------------------------------------------------*/
	
	public Organism fetchOrganism(int tid) {
		return index.get(tid);
	}
	
	public Sequence fetchSequence(int tid) {
		Organism o=index.get(tid);
		return o==null ? null : o.seqs.get(0);
	}
	
	public Organism fetchOrganism(String name) {
		TaxNode tn=tree.getNodeByName(name);
		return tn==null ? null : fetchOrganism(tn.name);
	}
	
	public Sequence fetchSequence(String name) {
		TaxNode tn=tree.getNodeByName(name);
		return tn==null ? null : fetchSequence(tn.name);
	}
	
	public ArrayList<SeqIndexResult> fetchSequence(Sequence s, SketchMakerMini smm, DisplayParams params) {
		if(params==null) {params=new DisplayParams();}
		if(smm==null) {
			params.format=DisplayParams.FORMAT_JSON;
			params.taxLevel=TaxTree.GENUS;
			smm=new SketchMakerMini(tool, SketchObject.ONE_SKETCH, 0, 0, (byte)0);
		}
		Sketch query=s.toSketch(smm, null);
		SketchResults results=sketchIndex.getSketches(query, params);
		
		CompareBuffer buffer=new CompareBuffer(false);
		int minHits=3;
		int maxHits=-1;
		for(Sketch b : results.refSketchList) {
			Comparison c=compareOneToOne(query, b, buffer, minHits, 0.0f, 0.6f);
			if(c!=null) {
				maxHits=Math.max(maxHits, c.hits());
				results.list.add(c);
				minHits=Math.max(minHits, (c.hits()/2)+1);
			}
		}
		Collections.sort(results.list);
		
		ArrayList<SeqIndexResult> alignments=new ArrayList<SeqIndexResult>(results==null ? 0 : results.list.size());
		IDAligner ida=new QuantumAligner();
		int position=0;
		for(Comparison c : results.list) {
//			if(c.hits()<minHits || ())//TODO
			Sketch refSketch=c.b;
			boolean is16S=(refSketch.r16S()!=null);
			byte[] ssu=(is16S ? refSketch.r16S() : refSketch.r18S());
			if(ssu==null) {continue;}
			float ani=ida.align(s.bases, ssu);
			Sequence seq=new Sequence(ssu, refSketch.taxName(), is16S ? "16S" : "18S", refSketch.taxID, refSketch.taxID);
			SeqIndexResult sir=new SeqIndexResult(s, seq, ani);
			sir.c=c;
			alignments.addLast(sir);
			position++;
		}
		Collections.sort(alignments);
		return alignments;
	}
	
	private static Comparison compareOneToOne(final Sketch a, final Sketch b, CompareBuffer buffer,
			int minHits, float minWKID, float minANI){
		if(a==b){return null;}
		final int matches=a.countMatches(b, buffer, null, false, null, -1);
		assert(matches==buffer.hits());
		if(matches<minHits){return null;}

		{
			final float wkid=buffer.wkid();
			if(wkid<minWKID){return null;}

			if(minANI>0){
				final float ani=buffer.ani();
				if(ani<minANI){return null;}
			}
		}
		Comparison c=new Comparison(buffer, a, b);
		return c;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Clade             ----------------*/
	/*--------------------------------------------------------------*/
	
	public int addClades(Iterable<Clade> clades) {
		int added=0;
		for(Clade c : clades) {
			if(c.taxID>0) {
				Organism o=index.get(c.taxID);
				if(o!=null && o.clade==null) {
					o.clade=c;
					added++;
				}
			}
		}
		return added;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Sketch            ----------------*/
	/*--------------------------------------------------------------*/
	
	public static SketchIndex sketchSeqs(Collection<Organism> organisms, String type) {
		ArrayList<Sequence> seqs=new ArrayList<Sequence>();
		for(Organism o : organisms) {
			if(type==null) {
				seqs.addAll(o.seqs);
			}else if(type.equalsIgnoreCase("16S")) {
				if(o.r16s!=null) {seqs.add(o.r16s);}
			}else if(type.equalsIgnoreCase("16S")) {
				if(o.r18s!=null) {seqs.add(o.r18s);}
			}else {assert(false) : type;}
		}
		BinObject.sketchDensity=density;
		BinObject.sketchInBulk=false;
		BinSketcher.send=false;
		BinSketcher bs=new BinSketcher(Shared.threads(), 32);
		bs.sketch(seqs, parseTID);
		
		ArrayList<Sketch> sketches=new ArrayList<Sketch>();
		for(Sequence seq : seqs) {
			if(seq.sketch!=null) {
				seq.sketch.taxID=seq.taxid();
				seq.sketch.name0=Sketch.fix(seq.name);
				if(seq.sketch.taxID>0) {
					TaxNode tn=tree.getNode(seq.tid);
					if(tn!=null) {seq.sketch.taxName=Sketch.fix(tn.name);}
				}
				sketches.add(seq.sketch);
				if(seq.r16S()) {seq.sketch.set16S(seq.bases);}
				else if(seq.r18S()) {seq.sketch.set18S(seq.bases);}
			}
		}
		
		SketchIndex sketchIndex=new SketchIndex(sketches);
		sketchIndex.load();
		return sketchIndex;
	}
	
	private static final SketchTool makeTool() {
		SketchObject.AUTOSIZE_LINEAR_DENSITY=density;
		SketchObject.AUTOSIZE_LINEAR=true;
		SketchObject.AUTOSIZE=false;
		SketchObject.SET_AUTOSIZE=true;
		SketchObject.minSketchSize=3;
		
//		SketchObject.AUTOSIZE=false;
//		SketchObject.defaultParams.minKeyOccuranceCount=2;
		SketchObject.defaultParams.parse("trackcounts", "trackcounts", null);
//		SketchObject.defaultParams.minProb=0;
		SketchObject.postParse();
		SketchObject.defaultParams.maxRecords=10;
		SketchObject.defaultParams.taxLevel=TaxTree.SPECIES;
		return new SketchTool(200, SketchObject.defaultParams);
		
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	public long readsProcessed;
	public long basesProcessed;
	public long sequencesAdded;
	public long organismsAdded;
	public long basesAdded;
	public HashMap<Integer, Organism> index;
//	public SketchIndex sketchIndex16S;
//	public SketchIndex sketchIndex18S;
	public SketchIndex sketchIndex;
	private static float density=0.125f;
	private final SketchTool tool=makeTool();

	/*--------------------------------------------------------------*/
	/*----------------           Statics            ----------------*/
	/*--------------------------------------------------------------*/
	
	public static boolean parseTID=true;
	public static int maxSeqsPerOrganism=Shared.MAX_ARRAY_LEN;
	
}
