-- This file works with an SQLite database generated from the contents
-- of backbone-current-simple.txt.gz

-- This file conflates the accepted usage key for synonyms with the parent key.
-- The query below shows that if the field is a synonym (or misapplied) then the
-- the parent_key always maps onto an accepted or doubtful taxon.

SELECT syn.status, par.status, count(*)
  FROM (SELECT status, parent_key
    FROM backbone
    WHERE status NOT IN ('DOUBTFUL','ACCEPTED')) AS syn
  INNER JOIN (SELECT id, status FROM backbone) as par
  ON syn.parent_key = par.id
  GROUP BY syn.status, par.status;

-- what counts as is_synonym?
SELECT
    status,
    count(CASE WHEN is_synonym = 't' THEN 1 END) as synonym,
    count(CASE WHEN is_synonym = 'f' THEN 1 END) as not_synonym
FROM backbone
GROUP BY status;


-- DOUBTFUL species

-- Have a parent in all but the pseudo-kingdom incertae sedis
select *
    from backbone
    where status = 'DOUBTFUL'
    and parent_key is null;

-- Find accepted species where the genus is doubtful: an example of higher
-- taxa not necessarily being accepted. Looking in Formicidae (4342)

select *
    from backbone as sp
    inner join (select * from backbone
                    where family_key = 4342
                    and rank = 'GENUS'
                    and status = 'DOUBTFUL') as dg
    on sp.parent_key = dg.id
    where sp.status = 'ACCEPTED'
    limit 10;


-- And one example: Liodes

select id, scientific_name, status, genus_key, family_key,
       order_key, class_key, phylum_key, kingdom_key
    from backbone
    where rank = 'GENUS'
    and canonical_name = 'Liodes';

-- A couple of different doubtful versions of it
select canonical_name, parent_key, family_key, order_key
    from backbone
    where parent_key in (3257142, 9016831)
    and rank = 'SPECIES'
    and status = 'ACCEPTED';

-- How common is duplicated scientific name (so no discriminating authority)
select scientific_name, count(*) as n
    from backbone
    where rank = 'SPECIES'
    and phylum_key = 54
    group by scientific_name
    order by n desc
    limit 10;

-- Example taxon
select * from backbone
    where canonical_name = 'Cancricepon choprae';

-- INCERTAE SEDIS
-- Taxa can have more than one rank and map down to kingdom Incertae Sedis

select distinct phylum_key, class_key, order_key, family_key, genus_key
    from backbone
    where kingdom_key = 0
    limit 10;

-- ... but none below family level
select count(distinct order_key)
    from backbone
    where kingdom_key = 0
    and order_key is not null;

select count(distinct family_key)
    from backbone
    where kingdom_key = 0
    and family_key is not null;

-- A handful of species have genus and family before dropping down to IS.
select count(*)
    from backbone
    where kingdom_key = 0
    and family_key is not null
    and genus_key is not null
    and rank = 'SPECIES';

-- OTHER QUERIES USED

-- Find species (ants) with non-genus level parents
SELECT a.id, a.canonical_name, b.canonical_name, b.rank
    FROM (
        SELECT id, canonical_name, parent_key, rank
            FROM backbone
            WHERE status == 'ACCEPTED') a
    JOIN (
        SELECT id, canonical_name, rank
            FROM backbone
            WHERE status == 'ACCEPTED') b
    ON a.parent_key = b.id
	WHERE a.rank == 'SPECIES'
	AND b.rank != 'GENUS'
	AND b.canonical_name == 'Formicidae'
	LIMIT 10;


select rank, count(*) from backbone
    where id in
        (select parent_key from backbone
            where rank = 'subspecies'
            and status != 'accepted'
            and parent_key is not null)
    group by rank;


select rank, count(*) from backbone
    where id in
        (select parent_key from backbone
            where rank = 'species'
            and status != 'accepted'
            and parent_key is not null)
    group by rank;


-- getting examples of oddly hooked in taxa
select * from backbone
    where id in
        (select parent_key from backbone
            where rank = 'form'
            and status != 'accepted'
            and parent_key is not null)
    and rank = 'variety'
    limit 10;

--
select * from backbone
    where rank = 'species'
    and status != 'accepted'
    and parent_key = 6410348;

select * from (
  SELECT
    canonical_name,
    class_key,
    count(CASE WHEN rank = 'SPECIES' THEN 1 END) as sp,
    count(CASE WHEN rank = 'SUBSPECIES' THEN 1 END) as ssp
  FROM backbone
    WHERE phylum_key = 44
    GROUP BY canonical_name) as nsp
  WHERE nsp.sp > 0 AND nsp.ssp > 0;

-- looking for cases where users might want to specify a particular
-- taxon rather than the accepted taxon that GBIF provides

SELECT spp.canonical_name
  FROM (SELECT canonical_name, status, parent_key
    FROM backbone
    WHERE canonical_name IN (
      -- find birds species names with two usages
      SELECT canonical_name FROM (
        SELECT
          canonical_name,
          count(*) as n
        FROM backbone
          WHERE class_key = 212
          and rank = 'SPECIES'
          GROUP BY canonical_name) as nsp
        WHERE nsp.n = 2)) AS spp
  INNER JOIN (SELECT id, rank FROM backbone) ataonttas par
  ON spp.parent_key = par.id
  WHERE par.rank = 'SUBSPECIES';


SELECT * FROM (
  SELECT
    canonical_name,
    count(*) as n
  FROM backbone
    WHERE class_key = 212
    and rank = 'SPECIES'
    GROUP BY canonical_name) as nsp
  WHERE nsp.n = 2;