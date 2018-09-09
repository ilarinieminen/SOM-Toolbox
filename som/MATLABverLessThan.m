
function result = MATLABverLessThan(version)

  % verLessThan() not currently (September 2018) available on Octave

  try
    result = verLessThan('matlab', version);
  catch
    result = true;
  end

end
