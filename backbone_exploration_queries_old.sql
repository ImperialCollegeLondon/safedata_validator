-- This file works with an SQLite database generated from the contents
-- of backbone-current.zip/Taxon.tsv

-- DOUBTFUL species

-- Never get and accepted usage...
select distinct acceptedNameUsageID 
    from backbone 
    where taxonomicStatus = 'doubtful';

-- ...but have a parent in all but the pseudo-kingdom incertae sedis
select *
    from backbone
    where taxonomicStatus = 'doubtful'
    and parentNameUsageID = "";

-- Find accepted species where the genus is doubtful: an example of higher
-- taxa not necessarily being accepted.

select *
    from backbone as sp
    inner join (select * from backbone
                    where taxonRank = 'genus'
                    and taxonomicStatus = 'doubtful') as dg
    on sp.parentNameUsageID = dg.taxonID
    where sp.taxonomicStatus = 'accepted'
    limit 10;

select *
    from backbone as sp
    inner join (select * from backbone
                    where family='Formicidae'
                    and taxonRank = 'genus'
                    and taxonomicStatus = 'doubtful') as dg
    on sp.parentNameUsageID = dg.taxonID
    where sp.taxonomicStatus = 'accepted'
    limit 10;


-- And one example: Liodes
select taxonID, scientificName, taxonomicStatus, genus, family, "order", class, phylum, kingdom
    from backbone
    where taxonRank = 'genus'
    and canonicalName = 'Liodes';

-- A couple of different doubtful versions of it
select canonicalName, parentNameUsageID, family, "order"
    from backbone
    where parentNameUsageID in (3257142, 9016831)
    and taxonRank = 'species'
    and taxonomicStatus = 'accepted';

-- How common is duplicated scientific name (so no discriminating authority)
-- I think this is where
select scientificName, count(*) as n
    from backbone
    where family = 'Formicidae'
    and taxonRank = 'species'
    group by scientificName
    order by n desc
    limit 10;

-- Example taxon
select * from backbone
    where canonicalName = 'Monomorium talbotae';

-- INCERTAE SEDIS
-- Taxa can have more than one rank and map down to kingdom Incertae Sedis

select distinct phylum, class, "order", family, genus
    from backbone
    where kingdom='incertae sedis'
    limit 10;

-- ... but none below family level
select count(distinct "order")
    from backbone
    where kingdom='incertae sedis'
    and "order" != "";

select count(distinct family)
    from backbone
    where kingdom='incertae sedis'
    and family != "";

-- A handful (52) of species have genus and family before dropping down to IS.
select count(*)
    from backbone
    where kingdom='incertae sedis'
    and family != ""
    and genus != ""
    and taxonRank = 'species'
    limit 10;

-- OTHER QUERIES USED

-- Find species (ants) with non-genus level parents
SELECT a.taxonID, a.canonicalName, b.canonicalName, b.taxonRank
    FROM (
        SELECT taxonID, canonicalName, parentNameUsageID, taxonRank
            FROM backbone
            WHERE taxonomicStatus == 'accepted') a
    JOIN (
        SELECT taxonID, canonicalName, taxonRank
            FROM backbone
            WHERE taxonomicStatus == 'accepted') b
    ON a.parentNameUsageID = b.taxonID
	WHERE a.taxonRank == 'species'
	AND b.taxonRank != 'genus'
	AND b.canonicalName == 'Formicidae'
	LIMIT 10;



select taxonRank, count(*) from backbone
    where taxonID in
        (select parentNameUsageID from backbone
            where taxonRank = 'subspecies'
            and taxonomicStatus != 'accepted'
            and parentNameUsageID is not null)
    group by taxonRank;


select taxonRank, count(*) from backbone
    where taxonID in
        (select parentNameUsageID from backbone
            where taxonRank = 'species'
            and taxonomicStatus != 'accepted'
            and parentNameUsageID is not null)
    group by taxonRank;


-- getting examples of oddly hooked in taxa
select * from backbone
    where taxonID in
        (select parentNameUsageID from backbone
            where taxonRank = 'form'
            and taxonomicStatus != 'accepted'
            and parentNameUsageID is not null)
    and taxonRank = 'variety'
    limit 10;

--
select * from backbone
    where taxonRank = 'species'
    and taxonomicStatus != 'accepted'
    and parentNameUsageID = 6410348;