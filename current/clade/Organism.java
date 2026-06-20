package clade;

import java.util.ArrayList;
import java.util.Iterator;

import bin.Sketchable;
import json.JsonObject;
import sketch.Sketch;
import sketch.SketchMakerMini;
import stream.Read;

public class Organism implements Iterable<Sequence>, bin.Sketchable {
	
	public Organism(String name_, int tid_) {
		name=name_;
		tid=tid_;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Methods           ----------------*/
	/*--------------------------------------------------------------*/
	
	public int add(Sequence s) {
		seqs.add(s);
		size+=s.size();
		if(s.r16S() && r16s==null) {r16s=s;}
		if(s.r18S() && r18s==null) {r18s=s;}
		return 1;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------          Sketchable          ----------------*/
	/*--------------------------------------------------------------*/
	
	@Override
	public int compareTo(Sketchable o){
		if(taxid()!=o.taxid()) {return taxid()-o.taxid();}
		//[clade/Organism#001] FIXED - was (int)(size()-o.size()): size is a long genome sum that can exceed 2^31 (large eukaryote), so the cast overflows -> wrong order. Long.compare is overflow-safe. Latent (Organism is used only by the dead SeqIndex subsystem). taxid()-o.taxid() is safe (taxIDs are small positive ints).
		return Long.compare(size(), o.size());
	}

	@Override
	public void setFrom(JsonObject jo){assert(false);}

	@Override
	public Sketch toSketch(SketchMakerMini smm, Read r) {
		if(sketch!=null) {return sketch;}
		if(r==null) {r=new Read(null, null, name, id());}
		r.id=name;
		r.numericID=id();
		for(Sequence c : seqs) {
			r.bases=c.bases;
			smm.processReadNucleotide(r);
		}
		return sketch=smm.toSketch(0);
	}

	@Override
	public void setID(int id){assert(false);}

	@Override
	public int id(){return tid;}

	@Override
	public float gc(){assert(false);return 0;}//unsupported-op fence (Organism has no single GC); same for setFrom/setID. Organism is dead (SeqIndex-only) -- these never fire.

	@Override
	public long size(){return size;}

	@Override
	public int taxid(){return tid;}

	@Override
	public int numContigs(){return seqs.size();}

	@Override
	public long sketchedSize(){return sketch==null ? 0 : sketch.length();}

	@Override
	public void clearTax(){}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	@Override
	public Iterator<Sequence> iterator(){return seqs.iterator();}
	
	public String name;
	public int tid;
	public long size;
	public Sequence r16s;
	public Sequence r18s;
	
	public ArrayList<Sequence> seqs=new ArrayList<Sequence>(2);
	public Sketch sketch;
	public Clade clade;
	
}
