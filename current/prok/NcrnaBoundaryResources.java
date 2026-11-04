package prok;

import ml.CellNet;

/** Loads the shared 10-input boundary adapter without altering its feature geometry.
 * @author Brian, Raiden */
public final class NcrnaBoundaryResources{

	/** Loads the selected family's network and both endpoint tables using its validated library. */
	public NcrnaBoundaryResources(final String family, final byte[][] library, final FamilyResourceResolver.Overrides overrides,
			final FamilyResourceResolver resolver){
		assert(library!=null && library.length>0 && resolver!=null)
			: "Shared boundary length features require the already validated family consensus library";
		// Every enabled boundary adapter requires all three roles, including 6S.
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
		meanLen=PackagedBoundaryConfig.upperMedianLength(library);
	}
	public final CellNet net;
	public final TrnaNinemerTableBuilder.LoadedTable start, stop;
	public final float meanLen;
}
