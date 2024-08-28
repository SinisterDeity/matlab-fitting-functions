function alphas = invertedOffsetLorentz(abc,angles)
alphas = (abc(1)-(abc(2).*abc(3))./((cosd(angles).^2)+((abc(3).^2).*sind(angles).^2)));
end