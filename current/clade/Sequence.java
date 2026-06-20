package clade;

import bin.Sketchable;
import json.JsonObject;
import sketch.Sketch;
import sketch.SketchMakerMini;
import stream.Read;
import structures.ByteBuilder;

public class Sequence implements bin.Sketchable {

	public Sequence(byte[] bases_, String name_, String type_, int tid_, long id_) {
		bases=bases_;
		name=name_;
		type=type_;
		tid=tid_;
		id=id_;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Methods           ----------------*/
	/*--------------------------------------------------------------*/
	
	public String toString() {
		return toBytes().toString();
	}
	
	public ByteBuilder toBytes() {
		ByteBuilder bb=new ByteBuilder();
		bb.append('>');
		if(name!=null) {bb.append(name);}
		bb.nl();
		//TODO: Possible bug [clade/Sequence#001] (LOW latent; dead code) - FASTA-wrap logic is broken/inverted: the wrap>0&&wrap<len branch (where wrapping SHOULD happen) appends the whole sequence unwrapped, while the else branch (incl. default wrap=-1) runs `for(i<bases.length) append(bases, i*wrap, wrap)` -> with wrap=-1 the offset i*wrap is negative -> AIOOBE; even with wrap>=len it runs len iterations reading past the end. The two branch bodies look swapped AND the loop should step i+=wrap (not i++) bounded to len. DEFERRED (a rewrite of an unused feature, not a typo): wrap is -1 by default with no caller, and toBytes()/toString() have no live callers (the whole SeqIndex subsystem is dead). Likely intended: if(wrap<=0||wrap>=len) append whole; else for(i=0;i<len;i+=wrap) append(bases,i,Math.min(wrap,len-i)).
		if(bases==null) {bb.nl();}
		else if(wrap>0 && wrap<bases.length) {
			bb.append(bases).nl();
		}else {
			for(int i=0; i<bases.length; i++) {
				bb.append(bases, i*wrap, wrap).nl();
			}
		}
		return bb;
	}

	public boolean r16S() {return "16S".equals(type);}
	public boolean r18S() {return "18S".equals(type);}
	public boolean ssu() {return r16S() || r18S() || "SSU".equals(type);}

	/*--------------------------------------------------------------*/
	/*----------------          Sketchable          ----------------*/
	/*--------------------------------------------------------------*/
	
	@Override
	public int compareTo(Sketchable o){
		if(taxid()!=o.taxid()) {return taxid()-o.taxid();}
		//(int)(size()-o.size()) is SAFE here (unlike Organism#001): size()=bases.length is int-bounded, so the long diff fits int. Dead anyway (SeqIndex-only).
		return (int)(size()-o.size());
	}

	@Override
	public void setFrom(JsonObject jo){assert(false);}
	
	@Override
	public Sketch toSketch(SketchMakerMini smm, Read r) {
		if(sketch!=null) {return sketch;}
		if(r==null) {r=new Read(null, null, name, id());}
		r.id=name;
		r.numericID=id();
		r.bases=bases;
		smm.processReadNucleotide(r);
		return sketch=smm.toSketch(0);
	}

	@Override
	public void setID(int id){assert(false);}

	@Override
	public int id(){return 0;}

	@Override
	public float gc(){assert(false);return 0;}

	@Override
	public long size(){
	return bases.length;}

	@Override
	public int taxid(){return -1;}

	@Override
	public int numContigs(){return 1;}

	@Override
	public long sketchedSize(){return sketch==null ? 0 : sketch.length();}

	@Override
	public void clearTax(){}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	public String name;
	public byte[] bases;
	public int tid;
	public final long id;
	
	public String type;
	public boolean amino=false;
	public Sketch sketch;

	/*--------------------------------------------------------------*/
	/*----------------            Statics           ----------------*/
	/*--------------------------------------------------------------*/
	
	public static int wrap=-1;
}
