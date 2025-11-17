package driver;

import java.util.ArrayList;
import java.util.Arrays;

import dna.AminoAcid;
import dna.Data;
import fileIO.ChainLine;
import shared.Shared;
import shared.Tools;
import var.VarLine;
import var.Variation;

/**
 * Translates genetic variation data between different genome builds.
 * Uses chain files to convert genomic coordinates and handles strand orientation.
 * Supports both VarLine and Variation object types for genomic variant data.
 * @author Brian Bushnell
 */
public class Translator {
	
	
	
	
	/**
	 * Creates a translator for converting between genome builds.
	 * Loads the appropriate chain files for coordinate translation.
	 * @param from_ Source genome build version
	 * @param to_ Target genome build version
	 */
	public Translator(int from_, int to_){
		fromBuild=from_;
		toBuild=to_;
		lines=Data.getChainLines(fromBuild, toBuild);
	}
	
	
	/**
	 * Translates an array of VarLine arrays between genome builds.
	 * Processes each VarLine individually and organizes results by chromosome.
	 * @param in Input array of VarLine arrays organized by chromosome
	 * @return Translated array of VarLine arrays, sorted by position within each chromosome
	 */
	public VarLine[][] translate(VarLine[][] in){
		ArrayList<VarLine>[] alvls=new ArrayList[in.length];
		for(int i=0; i<alvls.length; i++){
			alvls[i]=new ArrayList<VarLine>();
		}
		
		for(VarLine[] vla : in){
			if(vla!=null){
				for(VarLine vl : vla){
					VarLine vl2=translate(vl);
//					if(vl.haplotype==1 && (vl.intersects(244821744, 244821748) || (vl2!=null && vl2.intersects(246755120, 246755126)))){
//						System.out.println("\n"+vl+"\n->\n"+vl2);
//					}
					if(vl2!=null){
						int chrom=vl2.chromosome;
						alvls[chrom].add(vl2);
					}
				}
			}
		}
		
		VarLine[][] out=new VarLine[alvls.length][];
		for(int i=0; i<alvls.length; i++){
			out[i]=alvls[i].toArray(new VarLine[alvls[i].size()]);
			Arrays.sort(out[i]);
			alvls[i]=null;
		}
		
		return out;
	}
	
	
	/**
	 * Translates an array of Variation arrays between genome builds.
	 * Processes each Variation individually and organizes results by chromosome.
	 * @param in Input array of Variation arrays organized by chromosome
	 * @return Translated array of Variation arrays, sorted by position within each chromosome
	 */
	public Variation[][] translate(Variation[][] in){
		ArrayList<Variation>[] alvls=new ArrayList[in.length];
		for(int i=0; i<alvls.length; i++){
			alvls[i]=new ArrayList<Variation>();
		}
		
		for(Variation[] vla : in){
			if(vla!=null){
				for(Variation vl : vla){
					Variation vl2=translate(vl);
					if(vl2!=null){
						int chrom=vl2.chromosome;
						alvls[chrom].add(vl2);
					}
				}
			}
		}
		
		Variation[][] out=new Variation[alvls.length][];
		for(int i=0; i<alvls.length; i++){
			out[i]=alvls[i].toArray(new Variation[alvls[i].size()]);
			Arrays.sort(out[i]);
			alvls[i]=null;
		}
		
		return out;
	}
	
	
	/**
	 * Translates a single VarLine between genome builds.
	 * Uses chain files to convert coordinates and handles strand orientation.
	 * Applies reverse complement to sequence strings when translating to minus strand.
	 *
	 * @param v VarLine to translate
	 * @return Translated VarLine with updated coordinates and chromosome, or null if untranslatable
	 */
	public VarLine translate(VarLine v){
		
		ChainLine[] array=lines[v.chromosome];
		int index=ChainLine.binarySearch(v.beginLoc, array);
		if(index<0){return null;}
		ChainLine cl=array[index];
		if(!cl.contains(v.beginLoc, v.endLoc)){return null;}
		
//		System.out.println(cl);
		
		int[] dest1=cl.translate(v.beginLoc);
		int[] dest2=cl.translate(v.endLoc);
		
		if(dest1==null || dest2==null){return null;}
		
		VarLine v2=v.clone();
		
		assert(v!=null);
		assert(v2!=null) : v;
		
		if(cl.qStrand==Shared.PLUS){
			v2.chromosome=(byte)dest1[0];
			v2.beginLoc=dest1[2];
			v2.endLoc=dest2[2];
		}else{
//			assert(false) : "TODO";

			v2.chromosome=(byte)dest1[0];
			if(v.isPoint()){
				v2.beginLoc=v2.endLoc=dest1[2]-1;
			}else{
				v2.beginLoc=dest2[2];
				v2.endLoc=dest1[2];
			}
			
			if(v2.call!=null && Tools.isLetter(v2.call.charAt(0)) && !v2.call.equalsIgnoreCase("ref")){
				v2.call=AminoAcid.reverseComplementBases(v2.call);
			}
			
			if(v2.ref!=null && Tools.isLetter(v2.ref.charAt(0)) && !v2.ref.equalsIgnoreCase("ref")){
				v2.ref=AminoAcid.reverseComplementBases(v2.ref);
			}
			
		}
		
		assert(v2.endLoc-v2.beginLoc==v.endLoc-v.beginLoc) : "\n\n"+v.toSourceString()+"\n\n"+v2.toSourceString()+
		"\n\n"+v.beginLoc+" -> "+Arrays.toString(dest1)+
		"\n\n"+v.endLoc+" -> "+Arrays.toString(dest2)+
		"\n\n"+cl+"\n\n";
		
		assert(v2.beginLoc<=v2.endLoc) : "\n\n"+v.toSourceString()+"\n\n"+v2.toSourceString()+
		"\n\n"+v.beginLoc+" -> "+Arrays.toString(dest1)+
		"\n\n"+v.endLoc+" -> "+Arrays.toString(dest2)+
		"\n\n"+cl+"\n\n";
		
		v2.intern();
		return v2;
	}
	
	
	/**
	 * Translates a single Variation between genome builds.
	 * Delegates to VarLine translation if the object is a VarLine instance.
	 * Uses chain files to convert coordinates and handles strand orientation.
	 *
	 * @param v Variation to translate
	 * @return Translated Variation with updated coordinates and chromosome, or null if untranslatable
	 */
	public Variation translate(Variation v){
		
		if(v.getClass()==VarLine.class){
			return translate((VarLine)v);
		}
		assert(v.getClass()==Variation.class);
		
		ChainLine[] array=lines[v.chromosome];
		int index=ChainLine.binarySearch(v.beginLoc, array);
		if(index<0){return null;}
		ChainLine cl=array[index];
		if(!cl.contains(v.beginLoc, v.endLoc)){return null;}
		
		int[] dest1=cl.translate(v.beginLoc);
		int[] dest2=cl.translate(v.endLoc);
		if(dest1==null || dest2==null){return null;}
		
		Variation v2=v.clone();
		
		if(cl.qStrand==Shared.PLUS){
			v2.chromosome=(byte)dest1[0];
			v2.beginLoc=dest1[2];
			v2.endLoc=dest2[2];
		}else{
//			assert(false) : "TODO";

			v2.chromosome=(byte)dest1[0];
			if(v.isPoint()){
				v2.beginLoc=v2.endLoc=dest1[2]-1;
			}else{
				v2.beginLoc=dest2[2];
				v2.endLoc=dest1[2];
			}
			
			if(v2.call!=null && Tools.isLetter(v2.call.charAt(0)) && !v2.call.equalsIgnoreCase("ref")){
				v2.call=AminoAcid.reverseComplementBases(v2.call);
			}
			
			if(v2.ref!=null && Tools.isLetter(v2.ref.charAt(0)) && !v2.ref.equalsIgnoreCase("ref")){
				v2.ref=AminoAcid.reverseComplementBases(v2.ref);
			}
			
		}
		
		assert(v2.endLoc-v2.beginLoc==v.endLoc-v.beginLoc) : "\n\n"+v.toSourceString()+"\n\n"+v2.toSourceString()+
		"\n\n"+v.beginLoc+" -> "+Arrays.toString(dest1)+
		"\n\n"+v.endLoc+" -> "+Arrays.toString(dest2)+
		"\n\n"+cl+"\n\n";
		
		assert(v2.beginLoc<=v2.endLoc) : "\n\n"+v.toSourceString()+"\n\n"+v2.toSourceString()+
		"\n\n"+v.beginLoc+" -> "+Arrays.toString(dest1)+
		"\n\n"+v.endLoc+" -> "+Arrays.toString(dest2)+
		"\n\n"+cl+"\n\n";

		v2.intern();
		return v2;
	}
	
	
	/** Source genome build version */
	public final int fromBuild;
	/** Target genome build version */
	public final int toBuild;
	/** Chain file data for coordinate translation, organized by chromosome */
	public final ChainLine[][] lines;
	
}
