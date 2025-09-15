package var;

import java.io.File;
import java.util.ArrayList;

import align2.IndexMaker4;
import dna.ChromosomeArray;
import dna.Data;
import dna.FastaToChromArrays2;
import fileIO.ReadWrite;
import shared.Parse;
import shared.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import tracker.ReadStats;

/**
 * @author Brian Bushnell
 * @date Jul 23, 2012
 *
 */
public class ApplyVarsToReference {
	
	public static void main(String[] args){
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, new Object() { }.getClass().getEnclosingClass(), false);
			args=pp.args;
			//outstream=pp.outstream;
		}
		
		Timer t=new Timer();

		String inPattern=args[0];
		
		int minChrom=-1;
		int maxChrom=-1;
		int outgenome=-1;
		Data.GENOME_BUILD=-1;
		String name=null;
		
		for(int i=1; i<args.length; i++){
			final String arg=args[i].toLowerCase();
			String[] split=arg.split("=");
			String a=split[0];
			String b=split.length>1 ? split[1] : null;
			
			if(a.equals("ingenome")){
				Data.setGenome(Integer.parseInt(b));
				if(minChrom==-1){minChrom=1;}
				if(maxChrom==-1){maxChrom=Data.numChroms;}
			}else if(a.equals("outgenome")){
				outgenome=Integer.parseInt(b);
			}else if(a.equals("minchrom")){
				minChrom=Integer.parseInt(b);
			}else if(a.equals("maxchrom")){
				maxChrom=Integer.parseInt(b);
			}else if(a.equals("threads") || a.equals("t")){
				THREADS=Integer.parseInt(b);
			}else if(a.equals("nblocksize")){
				N_BLOCK_SIZE=Integer.parseInt(b);
			}else if(a.equals("nblocktrigger")){
				N_BLOCK_TRIGGER=Integer.parseInt(b);
			}else if(a.equals("staynearref")){
				STAY_NEAR_REF=Parse.parseBoolean(b);
			}else if(a.equals("append") || a.equals("app")){
				append=ReadStats.append=Parse.parseBoolean(b);
			}else if(a.equals("overwrite") || a.equals("ow")){
				overwrite=Parse.parseBoolean(b);
			}else if(a.startsWith("regen")){
				REGEN_N_BLOCKS=Parse.parseBoolean(b);
			}else if(a.startsWith("name=")){
				REGEN_N_BLOCKS=Parse.parseBoolean(b);
			}else{
				System.err.println("Unknown argument "+arg);
			}
		}
		
		assert(Data.GENOME_BUILD>-1);
		assert(outgenome>-1);
//		assert(Data.GENOME_BUILD!=outgenome);
		if(Data.GENOME_BUILD==outgenome){
			System.out.println("Warning! Overwriting input genome "+outgenome);
		}
		
		String fname=Data.chromFname(minChrom, outgenome);
		File f=new File(fname.substring(0, fname.lastIndexOf('/')));
//		assert(false) : f.getAbsolutePath();
		if(!f.exists()){f.mkdirs();}
		
		for(int chrom=minChrom; chrom<=maxChrom; chrom++){
			String outName=Data.chromFname(chrom, outgenome);
			assert(overwrite || !new File(outName).exists()) : "Destination "+outName+" already exists.";
//			assert(false) : inPattern+", "+outName;
			process(inPattern.replaceFirst("#", ""+chrom), outName, chrom);
		}
		
		FastaToChromArrays2.writeInfo(outgenome, maxChrom, (name==null ? Data.name : name), ""+Data.GENOME_BUILD+"_plus_variations", false, false);
		
		t.stop();
		
		{
			String path=IndexMaker4.fname(1, 1, 12, 1);
			int lastSlash=path.lastIndexOf('/');
			path=path.substring(0, lastSlash);
			File dir=new File(path);
			if(dir.exists()){
				System.out.println("Deleting old index for "+outgenome);
				for(File f2 : dir.listFiles()){
					if(f2.isFile() && (f2.getName().contains(".int2d") || f2.getName().endsWith(".txt"))){
						f2.delete();
					}
				}
			}
		}
		
//		System.out.println("Vars in: \t"+VARS_IN);
//		System.out.println("Vars out:\t"+VARS_OUT);
		System.out.println();
		System.out.println("Time: \t"+t);
		
	}
	
	public static void process(String inVarsName, String outChromName, int chrom) {
		ArrayList<Varlet> vars=Varlet.fromTextFile(inVarsName);
		ChromosomeArray cha=Data.getChromosome(chrom);
		ChromosomeArray chb=new ChromosomeArray(chrom, Shared.PLUS);
		
		//Next location to read in a
		int aloc=0;
		//Next location to set in b
		int bloc=0;
		
		for(int i=0; i<vars.size(); i++){
			
			Varlet v=vars.get(i);
			assert(v.beginLoc>=aloc) : i+"\n"+vars.get(i-1)+"\n"+v+"\n"; //Overlapping variations
			
			while(v.beginLoc<aloc){//skip it, for now.
				System.err.print("e");
				i++;
				if(i>=vars.size()){break;}
				v=vars.get(i);
			}
			
			if(STAY_NEAR_REF && Tools.absdif(aloc, bloc)>=REF_LIMIT){
				int dif=v.lengthDif();
				
				if(aloc<bloc){//skip insertions
					while(dif>0){
//						System.err.print("i");
						i++;
						if(i>=vars.size()){break;}
						v=vars.get(i);
						dif=v.lengthDif();
					}
				}else{//skip deletions
					while(dif<0){
//						System.err.print("d");
						i++;
						if(i>=vars.size()){break;}
						v=vars.get(i);
						dif=v.lengthDif();
					}
				}
			}
			
			//Advance to variation's beginning
			while(aloc<v.beginLoc){
				byte b=cha.get(aloc);
				chb.set(bloc, b);
				aloc++;
				bloc++;
			}
			
			//Apply variation
			if(v.varType==Variation.SNP){
				String call=v.call;
				String ref=v.ref;
				if(ref!=null && ref.equals("=")){ref=null;}
				for(int j=0; j<call.length(); j++){
					char c=call.charAt(j);
					if(ref!=null){
						assert(ref.charAt(j)==cha.get(aloc)) : "\n"+i+", "+v;
					}
					chb.set(bloc, c);
					aloc++;
					bloc++;
				}
			}else if(v.varType==Variation.DELINS){
				String call=v.call;
				for(int j=0; j<call.length(); j++){
					char c=call.charAt(j);
					chb.set(bloc, c);
					bloc++;
				}
				aloc+=v.lengthRef();
			}else if(v.varType==Variation.NOCALL){
				//Do nothing.  But, it should have been removed already.
				if(!foundNocall){
					System.err.println("*** Warning - found a nocall in input variations ***");
					foundNocall=true;
				}
			}else if(v.varType==Variation.NOREF){
				String call=v.call;
				for(int j=0; j<call.length(); j++){
					char c=call.charAt(j);
					assert(cha.get(aloc)=='N') : cha.get(aloc);
					chb.set(bloc, c);
					aloc++;
					bloc++;
				}
			}else if(v.varType==Variation.INS){
				String call=v.call;
				for(int j=0; j<call.length(); j++){
					char c=call.charAt(j);
					chb.set(bloc, c);
					bloc++;
				}
			}else if(v.varType==Variation.DEL){
				int len=v.lengthRef();
				assert(len>0);
				aloc+=len;
			}
		}
		
		//Finish writing array
		while(aloc<cha.array.length || aloc<=cha.maxIndex){
			byte c=cha.get(aloc);
			chb.set(bloc, c);
			aloc++;
			bloc++;
		}

		System.out.println("Length Shift for chr"+chrom+": \t"+(bloc-aloc));
		
		Data.unload(chrom, true);
		cha=null;
		
		if(REGEN_N_BLOCKS){
			chb=regenNBlocks(chb, N_BLOCK_SIZE, N_BLOCK_TRIGGER, N_BLOCK_END_SIZE);
		}
		
		chb.resize(chb.maxIndex+1);
		
		//Can't do this because it is read later
//		if(THREADS==1){ReadWrite.writeObjectInThread(cac, outChromName);}
//		else{ReadWrite.write(cac, outChromName);}
		
		ReadWrite.write(chb, outChromName, false);
	}
	
	public static ChromosomeArray regenNBlocks(ChromosomeArray cha, int blocksize, int trigger, int endsize){
		ChromosomeArray chb=new ChromosomeArray(cha.chromosome, cha.strand, cha.minIndex, cha.maxIndex);
		chb.maxIndex=-1;
		
		int aloc=0;
		int bloc=0;
		int ns=0;
		
		//Process start
		while(cha.get(aloc)=='N'){
			chb.set(bloc, 'N');
			ns++;
			aloc++;
			bloc++;
		}
		while(ns<endsize){
			chb.set(bloc, 'N');
			ns++;
			bloc++;
		}
		ns=0;
		
		
		//Process middle
		while(aloc<=cha.maxIndex){
			byte b=cha.get(aloc);
			if(b=='N'){
				ns++;
			}else{
				if(ns>=trigger){
					while(ns<blocksize){
						chb.set(bloc, 'N');
						bloc++;
						ns++;
					}
				}
				ns=0;
			}
			chb.set(bloc, b);
			aloc++;
			bloc++;
		}
		
		
		//Process end
		ns=0;
		for(int i=chb.maxIndex; i>=0; i--){
			if(chb.get(i)!='N'){break;}
		}
		while(ns<endsize){
			chb.set(chb.maxIndex+1, 'N');
			ns++;
		}
		
		return chb;
	}

	public static int THREADS=1;
	
	private static boolean foundNocall=false;
	private static boolean STAY_NEAR_REF=false;
	private static final int REF_LIMIT=20;
	public static boolean REGEN_N_BLOCKS=true;
	public static int N_BLOCK_END_SIZE=2000;
	public static int N_BLOCK_SIZE=300;
	public static int N_BLOCK_TRIGGER=80;
	/** Permission to overwrite existing files */
	public static boolean overwrite=true;
	/** Permission to append to existing files */
	public static boolean append=false;
	
}
