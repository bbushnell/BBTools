package prok;

import ml.CellNet;

/** Loads the shared 10-input boundary adapter without altering its feature geometry.
 * @author Brian, Raiden */
final class NcrnaBoundaryResources {

	NcrnaBoundaryResources(String family, byte[][] library, FamilyResourceResolver.Overrides overrides,
			FamilyResourceResolver resolver){
		assert(library!=null && library.length>0 && resolver!=null)
			: "Shared boundary length features require the already validated family consensus library";
		// Every enabled boundary adapter requires all three roles, including6S.
		// Family opt-outs are resolved before construction; missing selected resources fail loudly.
		final String netPath=resolver.resolve(family, "boundary", overrides).path;
		final String startPath=resolver.resolve(family, "boundary_start", overrides).path;
		final String stopPath=resolver.resolve(family, "boundary_stop", overrides).path;
		net=NcrnaBoundaryScorer.load(netPath);
		start=TrnaNinemerTableBuilder.loadTable(startPath);
		stop=TrnaNinemerTableBuilder.loadTable(stopPath);
		FamilyResourceManifest.require(start.type==TrnaBoundaryFeatures.BoundaryType.START,
			family+" boundary start table has type "+start.type+": "+startPath);
		FamilyResourceManifest.require(stop.type==TrnaBoundaryFeatures.BoundaryType.STOP,
			family+" boundary stop table has type "+stop.type+": "+stopPath);
		final int[] lengths=new int[library.length];
		for(int i=0; i<lengths.length; i++){lengths[i]=library[i].length;}
		java.util.Arrays.sort(lengths);
		meanLen=lengths[lengths.length/2];//CallGenes' historical upper-median feature.
	}
	final CellNet net;
	final TrnaNinemerTableBuilder.LoadedTable start, stop;
	final float meanLen;
}
