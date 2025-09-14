package stream;

import dna.Data;

/**
 * Transforms BBMap index coordinates into scaffold-relative coordinates.
 * @author Brian Bushnell
 * @date Aug 26, 2014
 *
 */
public class ScaffoldCoordinates {

	/*--------------------------------------------------------------*/
	/*----------------         Constructors         ----------------*/
	/*--------------------------------------------------------------*/
	
	public ScaffoldCoordinates(){}
	
	public ScaffoldCoordinates(Read r){set(r);}
	
	public ScaffoldCoordinates(SiteScore ss){set(ss);}
	
	/*--------------------------------------------------------------*/
	/*----------------            Methods           ----------------*/
	/*--------------------------------------------------------------*/
	
	public boolean set(Read r){
		valid=false;
		if(r.mapped()){setFromIndex(r.chrom, r.start, r.stop, r.strand(), r);}
		return valid;
	}
	
	public boolean set(SiteScore ss){
		return setFromIndex(ss.chrom, ss.start, ss.stop, ss.strand, ss);
	}
	
	public boolean setFromIndex(int iChrom_, int iStart_, int iStop_, int strand_, Object o){
		valid=false;
		if(iChrom_>=0){
			iChrom=iChrom_;
			iStart=iStart_;
			iStop=iStop_;
			if(Data.isSingleScaffold(iChrom, iStart, iStop)){
				assert(Data.scaffoldLocs!=null) : "\n\n"+o+"\n\n";
				scafIndex=Data.scaffoldIndex(iChrom, (iStart+iStop)/2);
				name=Data.scaffoldNames[iChrom][scafIndex];
				scafLength=Data.scaffoldLengths[iChrom][scafIndex];
				start=Data.scaffoldRelativeLoc(iChrom, iStart, scafIndex);
				stop=start-iStart+iStop;
				strand=(byte)strand_;
				valid=true;
			}
		}
		if(!valid){clear();}
		return valid;
	}
	
	public void clear(){
		valid=false;
		scafIndex=-1;
		iChrom=-1;
		iStart=-1;
		start=-1;
		iStop=-1;
		stop=-1;
		strand=-1;
		scafLength=0;
		name=null;
		valid=false;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields           ----------------*/
	/*--------------------------------------------------------------*/
	
	public int scafIndex=-1;
	public int iChrom=-1;
	public int iStart=-1, iStop=-1;
	public int start=-1, stop=-1;
	public byte strand=-1;
	public int scafLength=0;
	public byte[] name=null;
	public boolean valid=false;
	
}
